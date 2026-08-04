# tests/helpers.py
def strip_figure_text(fig):
    """
    Remove all text elements that vary across matplotlib versions,
    and lock figure geometry to prevent layout-engine-induced pixel
    shifts (e.g. ±1 px height differences between mpl versions).
    Use before returning a figure in mpl_image_compare tests.
    """
    # --- Remove all text artists ---
    fig.suptitle('')
    for text in fig.texts:
        text.set_text('')
    for ax in fig.axes:
        ax.set_title('')
        ax.set_xlabel('')      # not removed by remove_text=True
        ax.set_ylabel('')      # not removed by remove_text=True
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        leg = ax.get_legend()
        if leg:
            leg.remove()

    # --- Disable the layout engine ---
    # Prevents version-dependent reflow of subplot spacing at render/save time.
    try:
        fig.set_layout_engine("none")   # matplotlib >= 3.6
    except (AttributeError, ValueError):
        fig.set_tight_layout(False)
        try:
            fig.set_constrained_layout(False)
        except AttributeError:
            pass

    # --- Snap to integer pixel dimensions ---
    # Fractional inch*DPI values get rounded differently across versions,
    # producing ±1 px differences. Force exact integer pixel counts.
    dpi = fig.get_dpi()
    w_in, h_in = fig.get_size_inches()
    fig.set_size_inches(
        round(w_in * dpi) / dpi,
        round(h_in * dpi) / dpi,
    )

    # Force a redraw with locked geometry before pytest-mpl reads the buffer
    fig.canvas.draw()

    return fig