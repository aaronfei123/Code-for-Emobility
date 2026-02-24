# -*- coding: utf-8 -*-
"""
Voronoi-based mapping of truck flows to 132 kV buses (London) + visualization
— robust version with dedup, guard points, precision control, and fallbacks.

Inputs:
  - Substations CSV: columns must contain at least [id, lon, lat] in WGS84
  - Boundary: use BBOX or polygon file
  - 03_network-nodes.csv : columns [id, lon, lat] (or remap below)
  - 04_network-edges.csv : must contain endpoints (na, nb or remap), distances and flows
Outputs:
  - out/bus_energy_annual.csv
  - out/geo/london_voronoi_cells.geojson
  - out/geo/edge_segments.geojson (when USE_OVERLAY=True and there are intersections)
  - out/geo/edge_midpoints.geojson  (when USE_OVERLAY=False or no intersections)
  - out/plots/voronoi_subs_nodes.png
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon
from shapely.errors import TopologicalError
from scipy.spatial import Voronoi

# ---------- CONFIG ----------
SUB_CSV   = r"C:/2023Research/Emobility/transportationmap/london_132kv_substations_bigger2_main.csv"
NODES_CSV = r"C:/2023Research/Emobility/03_network-nodes.csv"
EDGES_CSV = r"C:/2023Research/Emobility/04_network-edges.csv"
BOUNDARY  = None  # polygon path (optional). If None, use BBOX.

# BBOX (south, west, north, east) — 如果没有边界文件，就用这个
BBOX = (51.28, -0.58, 51.75, 0.42)

OUT_DIR   = "out"
OUT_GEO   = os.path.join(OUT_DIR, "geo")
OUT_PLOT  = os.path.join(OUT_DIR, "plots", "voronoi_subs_nodes.png")
OUT_CSV   = os.path.join(OUT_DIR, "bus_energy_annual.csv")

# CRS
CRS_GEO    = "EPSG:4326"   # lon/lat
CRS_METRIC = "EPSG:27700"  # meters (GB National Grid)

# Energy params
KWH_PER_KM_URBAN = 1.05
KWH_PER_KM_RURAL = 0.95
KWH_PER_KM = (KWH_PER_KM_URBAN + KWH_PER_KM_RURAL) / 2.0
SEASON_MULT = 1.00

# Robustness params
SUBS_MERGE_TOL_M = 60.0   # 合并过近站点（米）
GUARD_PAD_M      = 3000.0 # 护城河宽度
PREC_GRID_M      = 0.1    # set_precision 网格
MICROBUF_M       = 0.01   # overlay 微缓冲

# Choose line→bus allocation method
USE_OVERLAY   = True      # True: overlay, False: sampling
SAMPLE_STEP_M = 300.0     # sampling 步长（米）

try:
    from shapely import set_precision as shapely_set_precision
except Exception:
    shapely_set_precision = None

def set_prec(g, grid=PREC_GRID_M):
    if shapely_set_precision:
        return shapely_set_precision(g, grid)
    return g  # fallback: do nothing

def ensure_dirs():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(OUT_GEO, exist_ok=True)
    os.makedirs(os.path.dirname(OUT_PLOT), exist_ok=True)

def read_boundary():
    if BOUNDARY and os.path.exists(BOUNDARY):
        bnd = gpd.read_file(BOUNDARY)
        if bnd.crs is None:
            bnd.set_crs(CRS_GEO, inplace=True)
        return bnd.to_crs(CRS_GEO)
    if BBOX:
        s, w, n, e = BBOX
        poly = Polygon([(w, s), (e, s), (e, n), (w, n), (w, s)])
        return gpd.GeoDataFrame(geometry=[poly], crs=CRS_GEO)
    raise FileNotFoundError("Provide BOUNDARY file path or set BBOX=(south, west, north, east).")

def merge_close_points_metric(gdf_metric, tol_m=SUBS_MERGE_TOL_M):
    pts = np.c_[gdf_metric.geometry.x, gdf_metric.geometry.y]
    used = np.zeros(len(pts), dtype=bool)
    keep_rows = []
    keep_xy = []
    for i in range(len(pts)):
        if used[i]:
            continue
        cluster = [i]
        used[i] = True
        for j in range(i+1, len(pts)):
            if not used[j] and np.linalg.norm(pts[j]-pts[i]) <= tol_m:
                used[j] = True
                cluster.append(j)
        xy = pts[cluster].mean(axis=0)
        keep_rows.append(gdf_metric.iloc[cluster[0]])
        keep_xy.append(xy)
    out = gpd.GeoDataFrame(keep_rows, crs=gdf_metric.crs).reset_index(drop=True)
    out["geometry"] = gpd.points_from_xy([p[0] for p in keep_xy], [p[1] for p in keep_xy], crs=out.crs)
    return out

def add_guard_points(bounds, pad=GUARD_PAD_M):
    minx, miny, maxx, maxy = bounds
    guards = np.array([
        [minx-pad, miny-pad], [maxx+pad, miny-pad], [maxx+pad, maxy+pad], [minx-pad, maxy+pad],
        [(minx+maxx)/2, miny-pad], [maxx+pad, (miny+maxy)/2],
        [(minx+maxx)/2, maxy+pad], [minx-pad, (miny+maxy)/2],
    ])
    return guards

def sample_line(line, step_m=300.0):
    if line is None or line.is_empty:
        return [], []
    L = float(line.length)
    if not np.isfinite(L) or L <= 0:
        return [], []
    n = max(1, int(np.ceil(L / step_m)))
    pts = [line.interpolate((i + 0.5) / n, normalized=True) for i in range(n)]
    seg_km = (L / n) / 1000.0
    return pts, [seg_km] * n


# ---------- Main ----------
def main():
    ensure_dirs()

    # 1) Read data
    subs = pd.read_csv(SUB_CSV)
    for a, b in [("Lon","lon"),("Longitude","lon"),("Lon.","lon"),
                 ("Lat","lat"),("Latitude","lat")]:
        if a in subs.columns and "lon" not in subs.columns and b=="lon":
            subs = subs.rename(columns={a:"lon"})
        if a in subs.columns and "lat" not in subs.columns and b=="lat":
            subs = subs.rename(columns={a:"lat"})
    assert {"id","lon","lat"}.issubset(subs.columns), "Substation CSV must have [id, lon, lat]."
    gdf_subs = gpd.GeoDataFrame(subs, geometry=gpd.points_from_xy(subs.lon, subs.lat), crs=CRS_GEO)

    bnd = read_boundary()
    # 1) 确保 id/lon/lat 正确
    nodes = pd.read_csv(NODES_CSV)
    cols = {c.lower(): c for c in nodes.columns}

    def pick(*ks):
        for k in ks:
            if k in cols: return cols[k]
        return None

    id_col = pick("id", "network_node_id", "node_id")
    lon_col = pick("lon", "longitude", "network_node_x", "x", "lon_dd")
    lat_col = pick("lat", "latitude", "network_node_y", "y", "lat_dd")
    assert id_col and lon_col and lat_col, "03 里找不到 id/lon/lat"
    nodes = nodes.rename(columns={id_col: "id", lon_col: "lon", lat_col: "lat"})
    nodes["id"] = pd.to_numeric(nodes["id"], errors="coerce").astype("Int64")

    # 2) 先按 BBOX 过滤“伦敦节点”
    south, west, north, east = BBOX
    nodes_ldn = nodes[(nodes.lat.between(south, north)) & (nodes.lon.between(west, east))].copy()
    if len(nodes_ldn) < 50:
        raise ValueError(f"伦敦 BBOX 内节点太少：{len(nodes_ldn)}")

    gdf_nodes = gpd.GeoDataFrame(nodes_ldn, geometry=gpd.points_from_xy(nodes_ldn.lon, nodes_ldn.lat), crs=CRS_GEO)

    # 读完 03 节点后（还未建边）
    print("NODES head:", nodes.head(3)[["id", "lon", "lat"]].to_dict("records"))
    print("NODES bounds (raw):", np.round([nodes["lon"].min(), nodes["lat"].min(),
                                           nodes["lon"].max(), nodes["lat"].max()], 4))

    # 3) generate edge
    edges = pd.read_csv(EDGES_CSV).rename(columns={
        "Network_Edge_ID": "edge_id",
        "Network_Node_A_ID": "na",
        "Network_Node_B_ID": "nb",
        "Distance": "dist_km",
        "Traffic_flow_trucks_2019": "trucks_2019",
        "Traffic_flow_trucks_2030": "trucks_2030",
    })
    for c in ["na", "nb", "edge_id"]:
        edges[c] = pd.to_numeric(edges[c], errors="coerce").astype("Int64")
    for c in ["dist_km", "trucks_2019", "trucks_2030"]:
        edges[c] = pd.to_numeric(edges[c], errors="coerce").fillna(0.0)

    keep_ids = set(gdf_nodes["id"].astype(int).tolist())
    edges = edges[edges["na"].isin(keep_ids) & edges["nb"].isin(keep_ids)].copy()

    nd = gdf_nodes.set_index("id")[["geometry"]]
    edges = edges.merge(nd, left_on="na", right_index=True, how="left") \
        .merge(nd, left_on="nb", right_index=True, how="left", suffixes=("_a", "_b"))

    def make_line(r):
        if pd.isna(r.geometry_a) or pd.isna(r.geometry_b): return None
        return LineString([r.geometry_a, r.geometry_b])

    edges["geometry"] = edges.apply(make_line, axis=1)
    gdf_edges = gpd.GeoDataFrame(edges.dropna(subset=["geometry"]), geometry="geometry", crs=CRS_GEO)

    print("Nodes in BBOX:", len(gdf_nodes))
    print("Edges kept (both ends in BBOX):", len(gdf_edges))

    print("EDGES bounds (pre-clip):", np.round(gdf_edges.total_bounds, 4))  # [minx, miny, maxx, maxy]

    # 2) Clip to boundary
    gdf_edges = gpd.clip(gdf_edges, bnd.buffer(0))

    print("Edges after clip:", len(gdf_edges))

    # 3) Build Voronoi (robust)
    subs_m = gdf_subs.to_crs(CRS_METRIC)
    subs_m = merge_close_points_metric(subs_m, tol_m=SUBS_MERGE_TOL_M)

    bnd_m  = bnd.to_crs(CRS_METRIC)
    bpoly  = bnd_m.unary_union.buffer(0)
    minx, miny, maxx, maxy = bpoly.bounds

    real_xy = np.array([(p.x, p.y) for p in subs_m.geometry])
    guards  = add_guard_points((minx, miny, maxx, maxy), pad=GUARD_PAD_M)
    vor_pts = np.vstack([real_xy, guards])

    # qhull options: more robust on collinear/degenerate configs
    vor = Voronoi(vor_pts, qhull_options="Qbb Qc Qx")

    # only keep regions for real points (exclude guard points)
    cells_polys = []
    cells_ids   = []
    for p_idx, region_idx in enumerate(vor.point_region[:len(real_xy)]):
        verts = vor.regions[region_idx]
        if not verts or -1 in verts:
            continue
        try:
            poly = Polygon(vor.vertices[verts]).buffer(0)
            poly = poly.intersection(bpoly)
            if not poly.is_empty:
                cells_polys.append(poly)
                cells_ids.append(int(gdf_subs.iloc[p_idx]["id"]))  # map back to original id order
        except Exception:
            continue

    gdf_cells_m = gpd.GeoDataFrame({"id": cells_ids},
                                   geometry=[set_prec(p, PREC_GRID_M) for p in cells_polys],
                                   crs=CRS_METRIC)

    # 4) Node→cell assignment (within → micro-coverage → nearest)
    nodes_m = gdf_nodes.to_crs(CRS_METRIC)
    cells4join = gdf_cells_m[["id","geometry"]].rename(columns={"id":"bus_id"})

    pt_join = gpd.sjoin(nodes_m, cells4join, how="left", predicate="within")
    null_mask = pt_join["bus_id"].isna()
    if null_mask.any():
        tiny = nodes_m.loc[null_mask].copy()
        tiny["geometry"] = tiny.buffer(0.05)  # 5 cm
        touch = gpd.sjoin(tiny, cells4join, how="left", predicate="intersects")
        pt_join.loc[touch.index.intersection(pt_join.index), "bus_id"] = touch["bus_id"].values
    null_mask = pt_join["bus_id"].isna()
    if null_mask.any():
        nn = gpd.sjoin_nearest(nodes_m.loc[null_mask], cells4join, how="left")
        pt_join.loc[nn.index, "bus_id"] = nn["bus_id"].values
    pt_out = pt_join.to_crs(CRS_GEO)[["id","bus_id","geometry"]]
    pt_out["lon"] = pt_out.geometry.x
    pt_out["lat"] = pt_out.geometry.y
    pt_out[["id","lon","lat","bus_id"]].sort_values(["bus_id","id"]).to_csv(
        os.path.join(OUT_DIR, "mapping_nodes_to_bus.csv"), index=False
    )

    # === 节点 to 母线 ===
    nodes_m = gdf_nodes.to_crs(CRS_METRIC)
    cells4join = gdf_cells_m[["id", "geometry"]].rename(columns={"id": "bus_id"})


    pt_join = gpd.sjoin(nodes_m, cells4join, how="left", predicate="within")


    null_mask = pt_join["bus_id"].isna()
    if null_mask.any():
        tiny = nodes_m.loc[null_mask].copy()
        tiny["geometry"] = tiny.buffer(0.10)  # 10 cm
        touch = gpd.sjoin(tiny, cells4join, how="left", predicate="intersects")
        pt_join.loc[touch.index.intersection(pt_join.index), "bus_id"] = touch["bus_id"].values


    null_mask = pt_join["bus_id"].isna()
    if null_mask.any():
        subs_nn = gdf_subs.to_crs(CRS_METRIC)[["id", "geometry"]].rename(columns={"id": "bus_id"})
        nearest_fix = gpd.sjoin_nearest(nodes_m.loc[null_mask], subs_nn, how="left")
        pt_join.loc[nearest_fix.index, "bus_id"] = nearest_fix["bus_id"].values

    # 导出
    pt_out = pt_join.to_crs(CRS_GEO)[["id", "bus_id", "geometry"]].copy()
    pt_out["lon"] = pt_out.geometry.x
    pt_out["lat"] = pt_out.geometry.y
    pt_out[["id", "lon", "lat", "bus_id"]].sort_values(["bus_id", "id"]).to_csv(
        os.path.join(OUT_DIR, "mapping_nodes_to_bus.csv"), index=False
    )
    print("Saved:", os.path.join(OUT_DIR, "mapping_nodes_to_bus.csv"))

    # 5) Edges to bus energy
    edges_m = gdf_edges.to_crs(CRS_METRIC)[["edge_id", "trucks_2019", "trucks_2030", "dist_km", "geometry"]].copy()
    print("edges_m rows after to_crs:", len(edges_m))

    # --- 1) clean---
    # 丢掉 None
    edges_m = edges_m[edges_m["geometry"].notna()].copy()

    edges_m = gpd.GeoDataFrame(edges_m, geometry="geometry", crs=CRS_METRIC)

    try:
        from shapely.validation import make_valid
        edges_m["geometry"] = edges_m.geometry.apply(lambda g: make_valid(g) if g is not None else None)
    except Exception:
        edges_m["geometry"] = edges_m.geometry.buffer(0)

    edges_m = edges_m[edges_m.geometry.notna()].copy()
    edges_m = edges_m[~edges_m.geometry.is_empty].copy()
    edges_m = edges_m[edges_m.geom_type.isin(["LineString", "MultiLineString"])].copy()

    edges_m = gpd.GeoDataFrame(edges_m, geometry="geometry", crs=CRS_METRIC)

    def _set_prec_safe(g):
        try:
            return set_prec(g, PREC_GRID_M)
        except Exception:
            return g

    edges_m["geometry"] = edges_m.geometry.apply(_set_prec_safe)

    # 最终长度统计
    edges_m["len_m"] = edges_m.length
    print("edge length stats (m):", edges_m["len_m"].describe())
    if edges_m["len_m"].count() == 0:
        raise RuntimeError("edge length 全 NaN")


    USE_OVERLAY = False
    hits = None

    # ---------- sampling  ----------
    def sample_line(line, step_m=SAMPLE_STEP_M):
        if (line is None) or line.is_empty:
            return [], []
        L = float(line.length)
        if not np.isfinite(L) or L <= 0:
            return [], []
        n = max(1, int(np.ceil(L / step_m)))
        try:
            pts = [line.interpolate((i + 0.5) / n, normalized=True) for i in range(n)]
        except Exception:
            pts = []
            for seg in getattr(line, "geoms", [line]):
                if seg.is_empty:
                    continue
                Ln = float(seg.length)
                if Ln <= 0:
                    continue
                k = max(1, int(np.ceil(Ln / step_m)))
                pts.extend([seg.interpolate((i + 0.5) / k, normalized=True) for i in range(k)])
            if not pts:
                return [], []
        seg_km = (L / n) / 1000.0
        return pts, [seg_km] * len(pts)

    samples = []
    for _, row in edges_m.iterrows():
        geom = row.geometry
        if (geom is None) or geom.is_empty:
            continue
        pts_list, seglens = sample_line(geom, step_m=SAMPLE_STEP_M)

        if not pts_list:
            try:
                mid = geom.interpolate(0.5, normalized=True)
                pts_list, seglens = [mid], [float(geom.length) / 1000.0]
            except Exception:
                continue
        for p, lk in zip(pts_list, seglens):
            samples.append({
                "edge_id": row.edge_id,
                "trucks_2019": row.trucks_2019,
                "trucks_2030": row.trucks_2030,
                "seg_len_km": lk,
                "geometry": p
            })

    print("sample points:", len(samples))
    if len(samples) == 0:

        mids = edges_m.geometry.interpolate(0.5, normalized=True)
        pts = gpd.GeoDataFrame(
            edges_m[["edge_id", "trucks_2019", "trucks_2030"]].assign(seg_len_km=edges_m.length / 1000.0),
            geometry=mids, crs=edges_m.crs
        )
    else:
        pts_df = pd.DataFrame(samples)
        pts = gpd.GeoDataFrame(
            pts_df.drop(columns=["geometry"]),
            geometry=gpd.GeoSeries(pts_df["geometry"], crs=edges_m.crs),
            crs=edges_m.crs
        )

    # 最近邻分配到变电站
    subs_for_nn = gdf_subs.to_crs(CRS_METRIC)[["id", "geometry"]].rename(columns={"id": "bus_id"})
    nn = gpd.sjoin_nearest(pts, subs_for_nn, how="left")
    hits = nn.rename(columns={"bus_id": "id"})  #

    # === 采样点映射导出为 CSV  ===
    hits.to_crs(CRS_GEO).to_file(os.path.join(OUT_GEO, "edge_midpoints.geojson"), driver="GeoJSON")
    hits.drop(columns="geometry", errors="ignore").to_csv(
        os.path.join(OUT_DIR, "mapping_edge_midpoints.csv"), index=False
    )

    # —— 采样结果去重压缩
    mid = hits.drop(columns="geometry", errors="ignore").copy()
    mid["points"] = 1
    mid_agg = (mid.groupby(["edge_id", "id"], as_index=False)
    .agg({
        "trucks_2019": "first",
        "trucks_2030": "first",
        "seg_len_km": "sum",
        "points": "sum"
    })
    )
    mid_agg.to_csv(os.path.join(OUT_DIR, "mapping_edge_to_bus_compact.csv"), index=False)
    print("Saved compact mapping:", os.path.join(OUT_DIR, "mapping_edge_to_bus_compact.csv"))

    print("Saved:", os.path.join(OUT_DIR, "mapping_edge_midpoints.csv"))

    # 导出
    try:
        hits.to_crs(CRS_GEO).to_file(os.path.join(OUT_GEO, "edge_midpoints.geojson"), driver="GeoJSON")
    except Exception as e:
        print("write edge_midpoints.geojson failed:", e)

    if "seg_len_km" not in hits.columns:
        hits["seg_len_km"] = hits.length / 1000.0

    for col in ["trucks_2019", "trucks_2030"]:
        hits[f"{col}_truck_km"] = pd.to_numeric(hits[col], errors="coerce").fillna(0.0) * hits["seg_len_km"]
        hits[f"{col}_MWh"] = hits[f"{col}_truck_km"] * KWH_PER_KM * SEASON_MULT / 1000.0

    if hits.empty or (hits["trucks_2019_MWh"].sum() == 0 and hits["trucks_2030_MWh"].sum() == 0):
        print("DEBUG first 5 pts:", pts.head() if len(pts) > 0 else "no pts")
        print("DEBUG subs extent (metric):", gdf_subs.to_crs(CRS_METRIC).total_bounds)
        raise RuntimeError("采样到的点没有分配到任何母线或 seg_len 全 0。请检查 pts 与 substations 是否同一坐标系/范围。")

    # 汇总
    agg_cols = {
        "trucks_2019_truck_km": "sum", "trucks_2030_truck_km": "sum",
        "trucks_2019_MWh": "sum", "trucks_2030_MWh": "sum"
    }
    bus_energy = hits.groupby("id").agg(agg_cols).reset_index()
    bus_energy["MWh_day_2019"] = bus_energy["trucks_2019_MWh"] / 365.0
    bus_energy["MWh_day_2030"] = bus_energy["trucks_2030_MWh"] / 365.0
    bus_energy.sort_values("id").to_csv(OUT_CSV, index=False)
    print("bus_energy rows:", len(bus_energy))

    # 6) Exports
    gdf_cells = gdf_cells_m.to_crs(CRS_GEO)
    gdf_cells.to_file(os.path.join(OUT_GEO, "london_voronoi_cells.geojson"), driver="GeoJSON")
    if USE_OVERLAY and "geometry" in hits.columns and hits.geom_type.isin(["LineString","MultiLineString","GeometryCollection"]).any():
        try:
            hits.to_crs(CRS_GEO).to_file(os.path.join(OUT_GEO, "edge_segments.geojson"), driver="GeoJSON")
        except Exception:
            pass

    tmp = pd.read_csv(os.path.join(OUT_DIR, "mapping_nodes_to_bus.csv"))
    print("mapping_nodes_to_bus rows:", len(tmp), " | bus_id 空值计数:", tmp["bus_id"].isna().sum())

    tmp2 = pd.read_csv(os.path.join(OUT_DIR, "mapping_edge_midpoints.csv"))
    print("mapping_edge_midpoints rows:", len(tmp2))

    # 7) Plot
    fig, ax = plt.subplots(figsize=(9, 9))
    bnd.boundary.plot(ax=ax, color="black", linewidth=1, label="Boundary")
    gdf_cells.plot(ax=ax, facecolor="lightgray", alpha=0.15, edgecolor="tab:gray", linewidth=1, label="Voronoi cell")
    gdf_subs.plot(ax=ax, color="red", markersize=30, marker="o",
                  edgecolor="white", linewidth=0.2, zorder=3, label="132 kV Substation")
    nodes_in = gpd.clip(gdf_nodes, bnd)
    nodes_in.plot(ax=ax, color="tab:blue", markersize=10, marker=".", zorder=2, label="Road network node")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    ax.set_title("Voronoi (132 kV) with Substations & Road Network Nodes")

    handles, labels = ax.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    ax.legend(uniq.values(), uniq.keys(), loc="upper right", frameon=True)

    minx, miny, maxx, maxy = bnd.total_bounds
    dx, dy = (maxx-minx)*0.005, (maxy-miny)*0.005
    ax.set_xlim(minx-dx, maxx+dx); ax.set_ylim(miny-dy, maxy+dy)

    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=220)
    plt.show()

    print("\nDone.")
    print("  - Energy table:", OUT_CSV)
    print("  - GeoJSON cells:", os.path.join(OUT_GEO, "london_voronoi_cells.geojson"))
    if USE_OVERLAY:
        print("  - Edge segments (if overlay used):", os.path.join(OUT_GEO, "edge_segments.geojson"))
    else:
        print("  - Edge midpoints (sampling):", os.path.join(OUT_GEO, "edge_midpoints.geojson"))
    print("  - Preview plot:", OUT_PLOT)

if __name__ == "__main__":
    main()
