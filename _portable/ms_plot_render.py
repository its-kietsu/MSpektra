"""Display-only spectrum rendering; numerical arrays are never changed."""
import numpy as np
from matplotlib.collections import LineCollection


class SpectrumSticks(LineCollection):
    """Cache a pixel envelope for raster views, retaining every stick for vectors.

    Dense positive sticks overlap from the same zero baseline. In each half
    pixel column the tallest actual sample preserves their visible envelope.
    Zoomed views with few samples draw every sample, without peak centroiding.
    """
    def __init__(self, x, y, **kwargs):
        x, y = np.asarray(x), np.asarray(y)
        valid = np.isfinite(x) & np.isfinite(y) & (y > 0)
        self._x, self._y = x[valid], y[valid]
        if len(self._x) > 1 and np.any(np.diff(self._x) < 0):
            order = np.argsort(self._x, kind='stable')
            self._x, self._y = self._x[order], self._y[order]
        self._view_key = None
        super().__init__([], **kwargs)

    def full_data(self):
        return self._x, self._y

    def _view_samples(self, raster):
        x, y = self._x, self._y
        if not raster or self.axes is None or not len(x):
            return x, y
        a, b = sorted(self.axes.get_xlim())
        start, end = np.searchsorted(x, [a, b], side='left')
        end = np.searchsorted(x, b, side='right')
        x, y = x[start:end], y[start:end]
        columns = max(100, int(np.ceil(self.axes.bbox.width * 2)))
        if len(x) <= 3 * columns or not b > a:
            return x, y
        if self.axes.get_xscale() == 'linear':
            edges = np.linspace(a, b, columns + 1)
            starts = np.unique(np.searchsorted(x, edges[:-1], side='left'))
            starts = starts[starts < len(x)]
        else:
            px = self.axes.transData.transform(np.column_stack((x, np.zeros(len(x)))))[:, 0]
            buckets = np.floor((px - px.min()) * 2).astype(np.int64)
            starts = np.r_[0, np.flatnonzero(np.diff(buckets)) + 1]
        # Keep the first true maximum per column, at its own measured m/z.
        run = np.repeat(np.arange(len(starts)), np.diff(np.r_[starts, len(x)]))
        maxima = np.maximum.reduceat(y, starts)
        picks = np.flatnonzero(y == maxima[run])
        picks = picks[np.r_[True, run[picks][1:] != run[picks][:-1]]]
        return x[picks], y[picks]

    def draw(self, renderer):
        if not self.get_visible():
            return
        raster = type(renderer).__name__ == 'RendererAgg'
        ax = self.axes
        key = (raster, tuple(ax.get_xlim()), ax.bbox.width, ax.get_xscale()) if ax else (raster,)
        if key != self._view_key:
            x, y = self._view_samples(raster)
            segments = np.zeros((len(x), 2, 2), dtype=float)
            segments[:, :, 0] = x[:, None]
            segments[:, 1, 1] = y
            self.set_segments(segments)
            self._view_key = key
        super().draw(renderer)


_CANVAS = None


def canvas_class():
    """Do not render invisible file pages during construction/DPI changes."""
    global _CANVAS
    if _CANVAS is None:
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        class VisibleCanvas(FigureCanvasWxAgg):
            def draw(self, drawDC=None):
                if not self.is_saving() and not self.IsShownOnScreen():
                    self._isDrawn = False
                    return
                super().draw(drawDC=drawDC)
        _CANVAS = VisibleCanvas
    return _CANVAS
