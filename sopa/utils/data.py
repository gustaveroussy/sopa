import logging

import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from shapely.geometry import Point, box
from spatialdata import SpatialData
from spatialdata.datasets import BlobsDataset
from spatialdata.models import Image2DModel, PointsModel, ShapesModel

log = logging.getLogger(__name__)


def uniform(
    *_,
    length: int = 2_048,
    cell_density: float = 0.01,
    n_points_per_cell: int = 50,
    n_genes: int = 5,
    c_coords: list[str] = ["DAPI", "CK", "CD3", "CD20"],
    sigma_factor: float = 0.4,
    seed: int = 0,
    save_vertices: bool = False,
    apply_blur: bool = True,
) -> SpatialData:
    """Generate a dummy dataset composed of cells generated uniformly in a square. It also has transcripts.

    Args:
        length: Size of the square, in pixels
        cell_density: Density of cells per pixel^2
        n_points_per_cell: Mean number of transcripts per cell
        n_genes: Number of gene names
        c_coords: Channel names
        sigma_factor: Factor used to determine `sigma` for the gaussian blur.

    Returns:
        A SpatialData object with a 2D image, the cells polygon boundaries, and the transcripts
    """
    np.random.seed(seed)

    grid_width = max(1, int(length * np.sqrt(cell_density)))
    dx = length / grid_width
    sigma = dx * sigma_factor
    n_cells = grid_width**2
    n_points = n_points_per_cell * n_cells

    log.info(
        f"Image of size ({len(c_coords), length, length}) with {n_cells} cells and {n_points_per_cell} transcripts per cell"
    )

    # Compute cell vertices (xy array)
    vertices_x = dx / 2 + np.arange(grid_width) * dx
    x, y = np.meshgrid(vertices_x, vertices_x)
    xy = np.stack([x.ravel(), y.ravel()], axis=1)
    xy += np.random.uniform(-dx / 2, dx / 2, size=xy.shape)
    xy = xy.clip(0, length - 1).astype(int)

    vertices = pd.DataFrame(xy, columns=["x", "y"])

    # Create image
    image = np.zeros((len(c_coords), length, length))
    image[0, xy[:, 1], xy[:, 0]] += 1
    if len(c_coords) > 1:
        image[np.random.randint(1, len(c_coords), len(xy)), xy[:, 1], xy[:, 0]] += 1
    if apply_blur:
        image = gaussian_filter(image, sigma=sigma, axes=(1, 2))
    image = (image / image.max() * 255).astype(np.uint8)

    # Create cell boundaries
    cells = [Point(vertex).buffer(sigma).simplify(tolerance=1) for vertex in xy]
    bbox = box(0, 0, length - 1, length - 1)
    cells = [cell.intersection(bbox) for cell in cells]
    gdf = gpd.GeoDataFrame(geometry=cells)

    # Create transcripts
    point_cell_index = np.random.randint(0, n_cells, n_points)
    points_coords = sigma / 2 * np.random.randn(n_points, 2) + xy[point_cell_index]
    points_coords = points_coords.clip(0, length - 1)
    df = pd.DataFrame(
        {
            "x": points_coords[:, 0],
            "y": points_coords[:, 1],
            "genes": np.random.choice([chr(97 + i) for i in range(n_genes)], size=n_points),
        }
    )
    points = {"transcripts": PointsModel.parse(df)}
    if save_vertices:
        points["vertices"] = PointsModel.parse(vertices)

    return SpatialData(
        images={"image": Image2DModel.parse(image, c_coords=c_coords, dims=["c", "y", "x"])},
        points=points,
        shapes={"cells": ShapesModel.parse(gdf)},
    )


def _to_mask(length: int, xy: list[tuple[int, int]], sigma: float):
    radius = int(sigma)
    circle_size = 2 * radius + 1
    circle = np.zeros((circle_size, circle_size), dtype=np.uint8)

    circle_y, circle_y = np.meshgrid(np.arange(circle_size), np.arange(circle_size))
    distance = np.sqrt((circle_y - radius) ** 2 + (circle_y - radius) ** 2)
    circle[distance <= radius] = 1

    mask = np.zeros((length, length))
    where = np.stack(np.where(circle), axis=1) - radius

    for i, (x, y) in enumerate(xy):
        where_i = where + np.array([y, x])
        where_i = where_i[
            (where_i[:, 0] >= 0)
            & (where_i[:, 1] >= 0)
            & (where_i[:, 0] < length)
            & (where_i[:, 1] < length)
        ].astype(int)
        mask[where_i[:, 0], where_i[:, 1]] = i + 1

    return mask


def blobs(
    *_,
    length: int = 1_024,
    n_points: int = 10_000,
    c_coords=["DAPI", "CK", "CD3", "CD20"],
    **kwargs,
) -> SpatialData:
    """Adapts the blobs dataset from SpatialData for sopa. Please refer to the SpatialData documentation"""
    _blobs = BlobsDataset(
        length=length, n_points=n_points, c_coords=c_coords, n_channels=len(c_coords), **kwargs
    )

    image = _blobs._image_blobs(
        _blobs.transformations,
        _blobs.length,
        _blobs.n_channels,
        _blobs.c_coords,
    )
    image.data = (image.data * 255).astype(np.uint8)

    points = _blobs._points_blobs(_blobs.transformations, _blobs.length, _blobs.n_points)
    genes = pd.Series(np.random.choice(list("abcdef"), size=len(points))).astype("category")
    points["genes"] = dd.from_pandas(genes, npartitions=points.npartitions)

    return SpatialData(images={"blob_image": image}, points={"blob_transcripts": points})
