from io import BytesIO

import matplotlib.pyplot as plt


def figure_to_bytes(
    figure,
    file_format: str,
    dpi: int = 600,
) -> bytes:
    """Convert a Matplotlib figure into downloadable bytes."""

    buffer = BytesIO()

    save_arguments = {
        "format": file_format,
        "bbox_inches": "tight",
        "facecolor": "white",
    }

    if file_format.lower() in {
        "png",
        "tiff",
        "tif",
        "jpeg",
        "jpg",
    }:
        save_arguments["dpi"] = dpi

    figure.savefig(
        buffer,
        **save_arguments,
    )

    buffer.seek(0)

    return buffer.getvalue()


def create_figure_exports(
    figure,
) -> dict:
    """Create publication-ready figure exports."""

    return {
        "preview_png": figure_to_bytes(
            figure,
            "png",
            dpi=180,
        ),
        "png_300": figure_to_bytes(
            figure,
            "png",
            dpi=300,
        ),
        "png_600": figure_to_bytes(
            figure,
            "png",
            dpi=600,
        ),
        "tiff_600": figure_to_bytes(
            figure,
            "tiff",
            dpi=600,
        ),
        "svg": figure_to_bytes(
            figure,
            "svg",
        ),
        "pdf": figure_to_bytes(
            figure,
            "pdf",
        ),
    }


def close_figure(
    figure,
) -> None:
    """Close a Matplotlib figure to release memory."""

    plt.close(figure)