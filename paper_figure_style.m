function paper_figure_style()
%PAPER_FIGURE_STYLE  Force print-ready graphics defaults for this repository.
%
%   The MATLAB desktop theme can default new figures to a dark background, which
%   prints badly and is unreadable in an IEEE two-column paper.  Every script
%   that generates a figure calls this once, immediately after `clear`, so the
%   generated .png/.pdf files are always white with black axes and text,
%   independent of the local MATLAB theme and font settings.
set(0, 'defaultfigurecolor', 'w');
set(0, 'defaultaxescolor',  'w');
set(0, 'defaultaxesxcolor', 'k');
set(0, 'defaultaxesycolor', 'k');
set(0, 'defaultaxeszcolor', 'k');
set(0, 'defaulttextcolor',  'k');
set(0, 'defaultlegendcolor', 'w');
set(0, 'defaultlegendtextcolor', 'k');
set(0, 'defaultaxesfontsize', 8);
set(0, 'defaulttextfontsize', 8);
set(0, 'defaultaxesgridalpha', 0.25);
set(0, 'defaultlinelinewidth', 1.2);
end
