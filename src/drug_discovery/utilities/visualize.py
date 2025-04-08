from deeporigin_molstar import JupyterViewer


def jupyter_visualization(func):
    def wrapper(*args, **kwargs):
        html_visualization = func(*args, **kwargs)
        return JupyterViewer.visualize(html_visualization)

    return wrapper
