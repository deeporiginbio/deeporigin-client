import os
import shutil
from itertools import islice

DEFAULT_MAX_DEPTH = 2


def move_file_with_extension(file_path, extension):
    dir_path = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    file_base_name = os.path.splitext(file_name)[0]
    target_file_path = os.path.join(dir_path, f"{file_base_name}.{extension}")

    if os.path.isfile(target_file_path):
        existing_files = [
            f
            for f in os.listdir(dir_path)
            if f.startswith(f"{file_base_name}_#") and f.endswith(f".{extension}")
        ]
        counter = len(existing_files) + 1
        new_file_path = os.path.join(
            dir_path, f"{file_base_name}_#{counter}.{extension}"
        )

        shutil.move(target_file_path, new_file_path)


def remove_file(file_path):
    if os.path.isfile(file_path):
        os.remove(file_path)


def chunker(iterable, size):
    iterator = iter(iterable)
    for first in iterator:
        yield [first] + list(islice(iterator, size - 1))


def calculate_box_min_max(box_center, box_dimensions):
    half_dims = [dim / 2 for dim in box_dimensions]
    min_corner = [center - half for center, half in zip(box_center, half_dims, strict=False)]
    max_corner = [center + half for center, half in zip(box_center, half_dims, strict=False)]
    return min_corner, max_corner


def calculate_box_dimensions(min_coords, max_coords):
    """
    Calculate box dimensions from minimum and maximum coordinates.

    Args:
        min_coords (List[float]): The minimum x, y, z coordinates of the box.
        max_coords (List[float]): The maximum x, y, z coordinates of the box.

    Returns:
        List[float]: The dimensions of the box [length, width, height].
    """
    # Ensure both inputs are lists or arrays of length 3
    if len(min_coords) != 3 or len(max_coords) != 3:
        raise ValueError("min_coords and max_coords must each have exactly 3 elements.")

    # Calculate the dimensions
    length = max_coords[0] - min_coords[0]
    width = max_coords[1] - min_coords[1]
    height = max_coords[2] - min_coords[2]

    return [float(length), float(width), float(height)]
