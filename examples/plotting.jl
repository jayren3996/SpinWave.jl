if !isdefined(@__MODULE__, :should_write_example_plots)

const PALETTE = ["#2864a6", "#c43c39", "#2f8f54", "#7b4ea3", "#d98c2b"]

should_write_example_plots() = get(ENV, "SPINWAVE_WRITE_EXAMPLE_PLOTS", "0") == "1"

function example_plot_dir()
    return get(
        ENV,
        "SPINWAVE_EXAMPLE_PLOT_DIR",
        normpath(joinpath(@__DIR__, "..", "docs", "src", "assets", "examples")),
    )
end

function write_dispersion_svg(filename, path, energies; title, ylabel="energy")
    x = collect(1:size(energies, 2))
    series = [vec(energies[mode, :]) for mode in axes(energies, 1)]
    labels = ["mode $mode" for mode in axes(energies, 1)]
    return write_series_svg(
        filename,
        x,
        series;
        labels=labels,
        title=title,
        ylabel=ylabel,
        xticks=path.ticks,
        xticklabels=path.labels,
        ymin=0.0,
    )
end

function write_series_svg(
    filename,
    x,
    series;
    labels,
    title,
    ylabel,
    xticks=[first(x), last(x)],
    xticklabels=["", ""],
    ymin=nothing,
    ymax=nothing,
)
    mkpath(dirname(filename))
    width = 760.0
    height = 430.0
    left = 78.0
    right = 34.0
    top = 54.0
    bottom = 68.0
    plot_w = width - left - right
    plot_h = height - top - bottom
    all_y = reduce(vcat, series)
    lo = ymin === nothing ? minimum(all_y) : ymin
    hi = ymax === nothing ? maximum(all_y) : ymax
    hi == lo && (hi = lo + 1)
    span = hi - lo
    lo = lo - 0.04span
    hi = hi + 0.08span

    xmin, xmax = extrema(x)
    xspan = xmax == xmin ? 1.0 : xmax - xmin
    yspan = hi - lo
    sx(v) = left + (v - xmin) / xspan * plot_w
    sy(v) = top + (hi - v) / yspan * plot_h

    open(filename, "w") do io
        _svg_header(io, width, height)
        println(io, """<rect width="100%" height="100%" fill="#ffffff"/>""")
        println(io, """<text x="$(width / 2)" y="28" text-anchor="middle" class="title">$(_escape(title))</text>""")
        _draw_axes(io, left, top, plot_w, plot_h)
        _draw_y_ticks(io, lo, hi, left, top, plot_w, plot_h)
        _draw_x_ticks(io, xticks, xticklabels, sx, left, top, plot_h)
        println(io, """<text x="$(width / 2)" y="$(height - 18)" text-anchor="middle" class="label">q-path</text>""")
        println(io, """<text x="22" y="$(top + plot_h / 2)" text-anchor="middle" transform="rotate(-90 22 $(top + plot_h / 2))" class="label">$(_escape(ylabel))</text>""")
        for (i, y) in pairs(series)
            points = join((_point(sx(x[j]), sy(y[j])) for j in eachindex(x)), " ")
            color = PALETTE[mod1(i, length(PALETTE))]
            println(io, """<polyline points="$points" fill="none" stroke="$color" stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round"/>""")
        end
        _draw_legend(io, labels, width - right - 126, top + 8)
        println(io, "</svg>")
    end
    return filename
end

function write_heatmap_svg(filename, path, grid; title, ylabel="energy")
    mkpath(dirname(filename))
    width = 760.0
    height = 430.0
    left = 78.0
    right = 34.0
    top = 54.0
    bottom = 68.0
    plot_w = width - left - right
    plot_h = height - top - bottom
    nω_full, nq_full = size(grid.intensity)
    q_indices = unique(round.(Int, range(1, nq_full; length=min(nq_full, 80))))
    ω_indices = unique(round.(Int, range(1, nω_full; length=min(nω_full, 60))))
    nω = length(ω_indices)
    nq = length(q_indices)
    maxval = maximum(grid.intensity)
    maxval > 0 || (maxval = 1.0)
    cell_w = plot_w / nq
    cell_h = plot_h / nω

    open(filename, "w") do io
        _svg_header(io, width, height)
        println(io, """<rect width="100%" height="100%" fill="#ffffff"/>""")
        println(io, """<text x="$(width / 2)" y="28" text-anchor="middle" class="title">$(_escape(title))</text>""")
        for (iq_out, iq) in pairs(q_indices), (iω_out, iω) in pairs(ω_indices)
            value = grid.intensity[iω, iq]
            t = sqrt(clamp(value / maxval, 0.0, 1.0))
            color = _heat_color(t)
            x = left + (iq_out - 1) * cell_w
            y = top + (nω - iω_out) * cell_h
            println(io, """<rect x="$(_fmt(x))" y="$(_fmt(y))" width="$(_fmt(cell_w + 0.2))" height="$(_fmt(cell_h + 0.2))" fill="$color"/>""")
        end
        _draw_axes(io, left, top, plot_w, plot_h)
        _draw_y_ticks(io, minimum(grid.omegas), maximum(grid.omegas), left, top, plot_w, plot_h)
        _draw_x_ticks(io, path.ticks, path.labels, v -> left + (v - 1) / max(nq_full - 1, 1) * plot_w, left, top, plot_h)
        println(io, """<text x="$(width / 2)" y="$(height - 18)" text-anchor="middle" class="label">q-path</text>""")
        println(io, """<text x="22" y="$(top + plot_h / 2)" text-anchor="middle" transform="rotate(-90 22 $(top + plot_h / 2))" class="label">$(_escape(ylabel))</text>""")
        println(io, "</svg>")
    end
    return filename
end

function _svg_header(io, width, height)
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(_fmt(width))" height="$(_fmt(height))" viewBox="0 0 $(_fmt(width)) $(_fmt(height))" role="img">""")
    println(io, """
<style>
text { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill: #222; }
.title { font-size: 20px; font-weight: 650; }
.label { font-size: 13px; }
.tick { font-size: 12px; fill: #444; }
.grid { stroke: #e5e7eb; stroke-width: 1; }
.axis { stroke: #333; stroke-width: 1.4; }
</style>""")
end

function _draw_axes(io, left, top, plot_w, plot_h)
    println(io, """<rect x="$(_fmt(left))" y="$(_fmt(top))" width="$(_fmt(plot_w))" height="$(_fmt(plot_h))" fill="none" class="axis"/>""")
end

function _draw_y_ticks(io, lo, hi, left, top, plot_w, plot_h)
    for t in range(0, 1; length=5)
        value = lo + t * (hi - lo)
        y = top + (1 - t) * plot_h
        println(io, """<line x1="$(_fmt(left))" x2="$(_fmt(left + plot_w))" y1="$(_fmt(y))" y2="$(_fmt(y))" class="grid"/>""")
        println(io, """<text x="$(_fmt(left - 10))" y="$(_fmt(y + 4))" text-anchor="end" class="tick">$(_fmt(value))</text>""")
    end
end

function _draw_x_ticks(io, ticks, ticklabels, sx, left, top, plot_h)
    for (tick, label) in zip(ticks, ticklabels)
        x = sx(tick)
        println(io, """<line x1="$(_fmt(x))" x2="$(_fmt(x))" y1="$(_fmt(top))" y2="$(_fmt(top + plot_h))" class="grid"/>""")
        println(io, """<text x="$(_fmt(x))" y="$(_fmt(top + plot_h + 24))" text-anchor="middle" class="tick">$(_escape(label))</text>""")
    end
end

function _draw_legend(io, labels, x, y)
    isempty(labels) && return
    legend_h = 22 * length(labels) + 10
    println(io, """<rect x="$(_fmt(x - 8))" y="$(_fmt(y - 16))" width="126" height="$(_fmt(legend_h))" rx="5" fill="#ffffff" stroke="#d1d5db"/>""")
    for (i, label) in pairs(labels)
        yy = y + (i - 1) * 22
        color = PALETTE[mod1(i, length(PALETTE))]
        println(io, """<line x1="$(_fmt(x))" x2="$(_fmt(x + 22))" y1="$(_fmt(yy))" y2="$(_fmt(yy))" stroke="$color" stroke-width="2.4"/>""")
        println(io, """<text x="$(_fmt(x + 30))" y="$(_fmt(yy + 4))" class="tick">$(_escape(label))</text>""")
    end
end

function _heat_color(t)
    lo = (247, 251, 255)
    mid = (107, 174, 214)
    hi = (8, 81, 156)
    if t < 0.5
        u = 2t
        return _rgb(_mix(lo, mid, u)...)
    end
    u = 2t - 1
    return _rgb(_mix(mid, hi, u)...)
end

_mix(a, b, t) = ntuple(i -> round(Int, (1 - t) * a[i] + t * b[i]), 3)
_rgb(r, g, b) = "#" * string(r, base=16, pad=2) * string(g, base=16, pad=2) * string(b, base=16, pad=2)
_fmt(x) = string(round(Float64(x); sigdigits=4))
_point(x, y) = "$(_fmt(x)),$(_fmt(y))"
_escape(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

end
