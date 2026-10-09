MSpektra 4.1 - portable Windows package   (09.10.2026)
======================================================

Analysis of LC-MS and HRMS data, with deconvolution:
  * LCMS Analysis   LC-MS and HPLC files (MS and PDA): Shimadzu .lcd,
                   Agilent .D, Waters .raw, Thermo .raw, mzML, mzXML, ANDI
  * HRMS Analysis   Bruker .d folders (maXis and other QTOF instruments:
                   analysis.baf or the newer analysis.tsf; timsTOF:
                   analysis.tdf), Agilent MassHunter .D folders, Waters
                   .raw folders, Thermo .raw files and mzML:
                   internal calibration, exact mass, formula finder
  * Deconvolute    the full deconvolution window for any mass spectrum

The deconvolution engine is UniDec 8.2.1 by Michael T. Marty
(https://github.com/michaelmarty/UniDec). LCMS Analysis, HRMS Analysis and
the Deconvolute window of the analysis windows run UniDec's engine built into
the library msengine.dll (compiled from UniDec's C source, BSD license, see
LICENSES\UniDec_LICENSE.txt); the UniDec windows of the start screen use the
unmodified UniDec 8.2.1 package with its own unidec.exe.

  >> If you publish results obtained with the deconvolution (in any of the
  >> three windows), please cite:
  >> Marty et al., Anal. Chem. 2015, 87 (8), 4370-4376.
  >> DOI: 10.1021/acs.analchem.5b00140

This folder contains the unmodified UniDec 8.2.1 code (from PyPI) running on
its own private copy of Python 3.12.10 (official python.org embeddable
build), plus the MSpektra add-ons in _portable. Nothing is installed on
the computer, no admin rights are needed, and it does not touch or depend on
Anaconda or any other Python on the PC.


VERSION HISTORY
---------------
The version is shown on the start screen, in the help window (?) of LCMS
and HRMS Analysis, and in the first line of each log file. It went up by
0.1 with every update until 3.3. Since 3.31 a small change adds 0.01, a
large change 0.1 and a big change 1.0.

4.1 (09.10.2026)
Files of every vendor open in their own raw formats.
* LCMS Analysis opens Agilent ChemStation/OpenLab .D and .dx, Waters .raw,
  Thermo .raw, mzML, mzXML and ANDI/AIA .cdf files (MS, PDA/UV or both,
  also HPLC files without MS).
* HRMS Analysis opens Agilent MassHunter .D (TOF, Q-TOF), Waters .raw and
  Thermo .raw (Orbitrap) files, besides Bruker .d and mzML.
* Sciex .wiff files give a short note to convert them to mzML with
  ProteoWizard MSConvert.
* Files dropped onto MSpektra.exe open in the window that fits their format.
* Reports name the instrument of every vendor.
* Agilent files are read with rainbow 1.5.3 (LGPL, see LICENSES), Waters
  files with rainbow (LCMS) or the Waters MassLynx library (HRMS), Thermo
  files with Thermo's RawFileReader; averaging, calibration and everything
  else run in the C++ core as for Shimadzu and Bruker data.
* Compare view figures in journal style: Arial, open frame (time axis and
  scale bar for stacked and offset traces), whole number time ticks,
  colour-blind safe colours, retention times clear of the other traces.

4.0 (09.10.2026)
* MS Analysis is now called MSpektra: MSpektra.exe, MSpektra.bat and
  MSpektra (console).bat start it. Saved analyses, method files and
  settings of earlier versions open as before.

3.77 (09.10.2026)
* LCMS Postrun and HRMS Postrun are now called LCMS Analysis and HRMS
  Analysis (windows, start screen, help, reports and this manual).

3.76 (09.10.2026)
* Shortened texts in the top bar (file name, calibration and others) end
  with "…" again instead of three broken characters.
* The deconvolution log is written as UTF-8, so names with special
  characters no longer stop it from being saved.

3.75 (09.10.2026)
The windows follow the screen they are on (laptops, scaled displays, large
monitors).
* Plots: on a laptop the chromatograms, spectra and PDA map get most of the
  window height (a plot is never lower than 42 % of the view, up to 280
  points); large monitors look as before. The tile size no longer needs
  changing when moving between screens.
* On small screens the top bar and the tool row stay on one line (labels
  give way to icons, last ones first), and the side panel, file list and
  status field are narrower.
* Pop up windows (calibration, deconvolution, mass shifts, formula finder,
  Graph properties, polymer, kinetics, report and the others) always open
  inside the screen, sized to their content, with their buttons in reach.
* Dialogs are no longer twice as wide on displays scaled above 100 %.
* A window saved on a monitor that is no longer connected opens on a present
  display.

3.65 (09.10.2026)
Fewer descriptions: explanatory hint lines removed and tooltips, notes,
dialogs, help windows and reports shortened everywhere.
Saved analyses
* Saved only after an edit and only next to the files edited; no Resume
  analysis.cmd is written.
* Opening a data file that has a saved analysis asks whether to resume it.
* Opening a saved analysis no longer changes your settings; mass shifts it
  uses that are missing from your list are added.
* A saved analysis that cannot be opened no longer stops saving or closing
  the other files.
* Two windows no longer save over the same analysis; the program version is
  stored in the file; Calibrate automatically stays as you set it.
* Faster start; a damaged saved analysis given at start shows a message, and
  data files given with it open too.
Compare
* Region areas are calculated in the C++ engine.
* With a region set, the table keeps the main peaks and guide line area % and
  adds the region areas; CSV and report show both.
* X values may be separated by commas, semicolons or new lines everywhere; a
  bad X value keeps the areas and says which value is wrong.
* A run without usable data in the region shows no data instead of 0 %.
* The chosen 100 % reference run comes back when it is ticked again.
* Undo records only applied region values and names its steps correctly.
* Peak times are named RT (min) everywhere and use the run's own time.
LCMS and HRMS Analysis
* Ctrl+O opens the raw data dialog directly; one Open menu in HRMS Analysis.
* The peak table says No peaks when there are none.
* Help windows are shorter and readable again, without garbled engine names
  such as "the Bayesian engine engine".
* Deconvolution errors, progress and quality notes no longer name UniDec.
Reports
* The HRMS report shows a deconvoluted mass with the decimals of the window.
* Explanations, long captions and notes removed; data, tables and settings
  kept.
Libraries and files
* mskinetics.dll and mspolymer.dll rebuilt with the same compiler as the
  engine (results unchanged).
* The rare crash during garbage collection after reading a file came from the
  old engine build and is fixed by the engine of 3.55.
* Unused code removed (among it the unused Python maximum entropy).

3.55 (08.10.2026)
Fixes from a full review of the program.
LCMS Analysis
* Files with two scan events of the same polarity (scan and SIM): each event
  has its own name, so peaks, integration and exported spectra no longer mix
  them up.
* Choosing a new time range while a spectrum is still being averaged no
  longer leaves the view showing one range and the spectrum another.
* The UV spectrum follows the time handle of the PDA map while you drag it.
* A single time typed in the Average fields shows the spectrum at that time;
  a background with only one time is reported instead of dropped silently.
* Integrating the same peak twice with Drag or Add range no longer halves its
  area %.
* Dragging the gap between two small tiles no longer shrinks the wrong tile.
HRMS Analysis
* Averaging profile scans whose m/z axes differ (zero trimmed profiles) no
  longer doubles peak heights or shifts the apex at high m/z.
* Averaged centroid spectra show one stick per ion, so the exact mass check
  reports the true error.
* The formula finder proposes perfluorinated and perchlorinated compounds such
  as PFUnDA and PFTeDA.
* The isotope match works for compounds with B, Fe, Se, Pt, Hg, Sn and similar
  elements, instead of always showing 0 %.
* With polarity switching, Calibrate on the selected range and Calibration
  details open the scans and calibrant of the right polarity.
* Bruker files are no longer locked after an open that failed or was
  cancelled, and an interrupted temporary copy of a write protected Bruker
  folder is no longer used as if it were complete.
* Fixed a crash or freeze when a large profile mzML file was read from two
  places at once (the engine library is now built on Windows).
Deconvolution and mass shifts
* Maximum entropy applies the minimum intensity once, so weaker species are no
  longer lost when the baseline is subtracted.
* The automatic peak width is correct on spectra with gaps of zeros (it could
  be about 40 times too wide).
* IsoDec works when the program is in a folder whose name has accented or
  other non English characters.
* Mass shifts: conjugates of isotope resolved proteins with heavy tags are
  recognised, so the degree of conjugation counts every species. The table
  shows the same isotope mass as the mass table, the plot label and the saved
  peak list. The degree of conjugation by area is left empty, with a note,
  for methods that give no peak areas.
* Every result keeps its own files (single scans are named by their time; a
  number is added instead of writing over an earlier result), and edits only
  rewrite a result's own peak list.
* The fit is drawn only on the spectrum it was made from, and Deconvolute
  again offers the result's own spectrum.
Compare
* In the Offset layout, Align and Scale to a peak use the peak you click, not
  a time moved by the skew.
* X values typed for the area table stay with their runs when the list is
  sorted, moved or reversed.
* The y axis, scale bar and exported traces say % only when the traces were
  really scaled.
* λmax uses the peak shown in the time window when the runs are aligned.
* A short gap with no values in a trace no longer spreads over the rest of the
  trace when the rolling baseline is used.
Reports and import
* Deconvolution reports show the masses with as many decimals as the window.
* Expected masses are shown as typed, and 148,056.3 is read as 148056.3.
* HRMS and SI reports say "outside tolerance" when the checked formula's
  nearest peak lies outside the tolerance.
* PDF reports no longer fail with very long notes, and long sample names are
  shortened in the page header.
* JCAMP-DX import: correct m/z axis for compressed data without FIRSTX and
  LASTX and on NTUPLES pages, point count checked against NPOINTS, and times
  in milliseconds read correctly.
Saving, undo and settings
* Opening a saved analysis no longer adds an undo step: Ctrl+Z right after
  opening no longer throws the restored results away.
* Opening a saved analysis no longer replaces your deconvolution and
  calibration settings.
* Settings are no longer lost or reset when several windows save them at the
  same time.
* Keys pressed in the Polymer and Kinetics windows act in those windows, not
  on the main window.
* Undo in one window no longer switches the link MS and PDA times choice made
  in another window.
* After a failed recovery, the last good checkpoint is kept as a separate file
  that later saves do not replace.
* The launcher's error messages name MSpektra and its real start files.

3.45 (08.10.2026)
* Compare: with the Offset layout or a skew, the selected region is
  integrated over the same retention time in every run. Before, each run
  moved right by the skew was integrated over a window moved left by the
  same amount, so relative areas and kinetics compared different parts of
  the runs. The fill of each run is drawn under its own peak, moved with
  its trace; the dashed lines mark the region on the time axis.

3.44 (08.10.2026)
* Compare: the 100 % area reference and the line width and style set for
  each trace in Graph properties are kept when a saved analysis is opened
  again. Before, the reference went back to "largest area (automatic)"
  without notice and the trace lines to their defaults.

3.43 (08.10.2026)
* Compare: Excel, CSV, the report and kinetic fitting use the region start
  and end as typed in the fields, also without pressing Enter. Before, they
  used the previous region until Enter was pressed.

3.42 (08.10.2026)
* The analysis is saved also when the raw data is on another drive or a
  network share than the analysis folder. Before, saving failed and the file
  and the window could not be closed.
* When the analysis cannot be saved (for example a folder without write
  access), closing a file or the window asks "Close without saving the
  analysis?" instead of refusing to close.

3.41 (08.10.2026)
* A saved analysis with a polymer analysis window open can be reopened
  again. Before, reopening it failed, autosave stopped and the files could
  not be closed.

3.4 (08.10.2026)
* Automatic analysis saving and reopening, with a Save button and
  Ctrl+Shift+S, and an Open menu (raw data, analysis folder, project).
* Compare: region areas in every run with a 100 % reference, area vs X plot,
  Excel export and kinetic fitting.
* HRMS Analysis: polymer analysis (Mn, Mw, dispersity, repeat units, chain
  lengths) and the isotope details tile.
* Maximum entropy in its own library (msmaxent.dll).
* Plain HRMS and deconvolution reports.
* Faster stick spectra; the toolbar wraps on narrow windows.
* Updated display names in the workspace and reports.

3.35 (04.10.2026)
* Mass shifts on deconvolution results (LCMS Analysis and HRMS Analysis):
  right click a zero charge mass spectrum > Mass shifts (or Mass shifts...
  in the Deconvolution section of the side panel). Brackets above the peaks
  join the masses that differ by a known shift (Na, K, oxidation, water,
  carbamylation from urea, TFA, hexose and others; the list can be edited)
  or by one or more tags (a dye or a PDI tag, set by you) plus up to two
  other shifts (e.g. tag + water: the hydrolysed succinimide), each labelled
  with its name and the difference, e.g. "PDI + water +861.0"; a mass whose
  difference from the reference matches nothing is marked above its peak,
  e.g. "+429.1 unknown". A table next to the plot lists the pairs (from, to,
  difference, match, error, what the mass is) and, with a tag, the degree of
  conjugation: the share of the species with 0, 1, 2, ... tags by height and
  by area and the tags per molecule. The reference is chosen by itself (the
  tallest mass, or with a tag the lightest of the tag series when a
  conjugate is the tallest); right click a peak > Use ... Da as the
  unmodified species changes it. Proteins are compared by average masses,
  isotope resolved results by monoisotopic differences. The brackets are
  part of copied and saved images, the tables part of the reports. Computed
  in the C++ library.
* Method presets in LCMS Analysis and HRMS Analysis: Method in the top bar
  shows the method of the file shown; its menu saves every setting of the
  window under a name, applies a method to the file shown or to every open
  file, updates, renames, deletes, exports and imports methods
  (name.msmethod.json) and sets the default method, applied to every file
  opened in that window. A method includes the settings of the mass
  shift finder (shift list, tags, tolerances). See METHOD PRESETS below.
* HRMS Analysis, formula finder: it opens with the element limits and the
  isotope ranking used last in that file (or set by a method) instead of
  the defaults every time, and with the ion chosen under Exact mass when
  that ion has the polarity of the peak (as Check a formula does).
* Undo and Redo in LCMS Analysis and HRMS Analysis: every edit of a window
  can be undone (Ctrl+Z) and redone (Ctrl+Y or Ctrl+Shift+Z), up to 100
  steps per window. Undo and Redo are also in the top bar (greyed when
  there is nothing to undo; the tooltip names the step, e.g. "Undo: delete
  peak 9.04 min") and at the end of the right click menu of every plot.
  Undo shows the file and the view where the step was made and says in
  the status bar what was undone. Applying a method is one step, and the
  mass shift finder (on or off, reference, settings) is covered too. See
  UNDO below for what is covered.
* Keyboard shortcuts everywhere in the analysis windows, from one list:
  F1 opens it (searchable, by window). Ctrl combinations work anywhere in
  the window (Ctrl+O, Ctrl+W, Ctrl+R as before; new: Ctrl+1/2/3 for the
  views, Ctrl+E export the data of the plot clicked, Ctrl+D deconvolute,
  Ctrl+I integrate, Ctrl+Page Down / Page Up for the next or previous
  file); single letters choose the tools while a plot has the focus (S
  Select, Z Zoom, P Pan, B Background, I, D, K, R integration tools, X Mass
  chrom., U Measure, L Label, F Formula (HRMS), A Align, G Guide line, M
  m/z (Compare)); Home full view, + and - zoom, the arrow keys move the
  picked scan (Mass spectrometry view) or the sliders (PDA view), Delete
  deletes the peak chosen. Typing in a field is never taken as a shortcut.
  The tooltips of the tool bar and the menus show the keys. See KEYBOARD
  SHORTCUTS below.

3.34 (03.10.2026)
* Compare view, m/z tool: the two spectrum tiles work like the spectra of
  the Mass spectrometry view: the m/z under the pointer, drag to zoom the
  m/z axis, Ctrl + drag or the Zoom tool for a box, Ctrl + wheel, Pan,
  double click for the full view. While m/z is on, the tool bar also
  offers the spectrum tools of the Mass spectrometry view: Mass chrom.
  (click a peak: every run shows the mass chromatogram of that ion),
  Measure (two peaks: their distance, the isotope spacing or, for
  adjacent charge states, the charges and the mass, for the polarity of
  the tile) and Label (pin the label of a peak). Right click a spectrum:
  mass chromatogram of that ion in every run, label it or remove every
  label, remove the measurement, display as profile or sticks, show it in
  the Mass spectrometry view, export it (text, JCAMP-DX or CSV with the
  run and times), copy the data, open it in the Deconvolute window, full
  view, copy or save the image.
* Compare view, m/z tool: no adduct names under the peaks any more; the
  labels are those of the Mass spectrometry view. The suggested mass in
  the line above the tiles is unchanged.
* Compare view, m/z tool: regions on the traces, as on the chromatograms
  of the Mass spectrometry view. While m/z is on, drag across a trace to
  average that time range of its run (both polarities), Shift + drag to
  set the background range of that run (subtracted from its spectra; for
  a peak clicked it replaces the spectrum just before the peak while the
  side panel option is on), Alt + click for one scan at that time. A
  click still picks the peak, Ctrl + drag zooms, Zoom and Pan work as
  before. The averaged range and the background range are shaded on that
  trace only (the background grey and hatched), and the tiles and the
  report say e.g. "averaged 9.20 to 9.60 min (range), background 8.80 to
  9.00 min". Right click the plot to remove a background range; m/z off
  removes them. With m/z off a drag zooms the time axis as before.
* The drawing, labels, measuring, export and copy of the spectra are one
  piece of code for the Mass spectrometry view and the Compare view (no
  change in the Mass spectrometry view).

3.33 (03.10.2026)
* Compare view: m/z in the tool bar (on or off) shows the mass spectra of
  a peak. Click a peak of any trace: the positive and negative ion spectra
  of that run, averaged over the scans of the peak (or one scan at its
  top), appear in two tiles below the plot, with the strongest ions
  labelled and a suggested neutral mass from the adducts found in both
  polarities (computed in the C++ library; mobile phase acids such as TFA
  named as such, Br and Cl M+2 patterns noted, unexplained base peaks
  named, a weak case called weak). PDA peaks are taken at the PDA time
  plus the MS detector delay of the file; aligned, shifted or skewed
  traces at the time of the run itself. The peak is marked on its trace.
  Right click a spectrum to follow that ion through every run (mass
  chromatogram), to open it in the Mass spectrometry view, or to export it
  as CSV. The Comparison report can include the two spectra with their
  caption. Spectra in the side panel: averaging, background before the
  peak, number of ions labelled.

3.32 (03.10.2026)
* Compare view, new layout Offset (2D waterfall): each run a little higher
  and to the right of the next one, overlapping, with a real y axis with
  its values, the time axis of the runs themselves, the front run drawn
  over the ones behind it and the name of each run right of the frame in
  its colour (the look of stacked HPLC traces in papers). Offset keeps its
  own spacing and skew (15 % and 4 % to start with).
* Compare view: a Y axis choice (automatic, scale bar, axis with values,
  none) instead of the Scale bar check box, and names "Right of the frame"
  for any layout. Long names there go on more lines in narrow images.
* Compare view: right click > Graph properties: font (e.g. Times New Roman),
  text sizes and bold, axis titles, x and y range, tick steps, ticks inside
  or outside, minor ticks, frame (box, left and bottom, bottom only, none),
  grid, colour, width and style (solid, dashed, dotted) of each trace, fill
  opacity, legend position and frame, a title above the plot and a panel
  letter (A, B). Apply shows the change, Cancel undoes it, Reset all goes
  back to the plot as drawn without them. They last for the session (not
  saved) and apply to copied and saved images and the Comparison report.
* Copied and saved images: labels are placed once more when they changed
  the height of the plot (a title on more lines than on screen).

3.31 (03.10.2026)
* Compare view: its calculations now run in the C++ library (msengine.dll)
  like the rest of the program: blank subtraction, the peak search of the
  alignment, offset and drift baselines, the processing of each trace
  (smoothing, shift, time window, baseline, scale), the stacking offsets,
  the numbers of the peak table (main peak, area %, area % at the guide
  lines) and lambda max of the main peak. The results are the same as
  before to the last digit (checked on 16,000 cases); the Python code stays
  as the fallback.
* New installation folder: MSpektra (next to MSpektra 1.8, which is
  left as it is). Unpacking it copies the settings of MSpektra 1.8
  (config\portable_settings.json); from then on each folder keeps its own.

3.3 (03.10.2026)
* LCMS Analysis has a third view, Compare: the chromatograms of every open
  file on top of each other (stacked, for presentations and reaction
  monitoring) or overlaid. See COMPARE VIEW below.
* Reports: a new Comparison report (PDF or Word) of the Compare view: the
  plot exactly as shown, the signal and every processing step, a table of
  the runs (sample, date, method, vial, injection volume, time shift), the
  peaks of each trace (main peak, height, area %, area % at each guide
  line), notes when the runs differ, and a methods text. The report window
  scrolls when it would be taller than the screen; a report that is not
  available shows the reason in one line.
* Word reports: the columns of the tables keep their widths in LibreOffice
  too (they came out equal).

3.2 (03.10.2026)
* Reports show what the window shows. Every figure that has a tile in the
  window (chromatograms, mass spectra, UV spectrum, zero charge mass spectra,
  charge states) is now the image of that tile, drawn for the page: the same
  view (zoom), traces, labels (also Label all peaks in view, pinned and
  added masses), markers, the fit and charge states on the spectrum, shaded
  time ranges and mass chromatograms. Before, the reports drew figures of
  their own from the data (e.g. the deconvolution over another mass range,
  with other labels, and the charge states even when hidden in the window).
  The trace names of a tile are written into the plot.
* HRMS report: every result shown in the window is in it, newest first, each
  with its notes and its mass table; the charge states only when that tile
  is shown; hidden results stay out. Without a checked formula, the exact
  mass table lists the peaks labelled on the spectrum in the window.
  Deconvolution report: the spectrum tile with the fit and the charge states
  of the mass selected in the table, and the result as its row shows it.
  Supporting Information page: its panels are the tiles too (the isotope
  pattern, which has no tile, is drawn as before).
* A page is narrower than a tile on a large screen: the text of such a
  figure is then up to 20 % smaller than elsewhere in the report, so that
  the labels placed for the tile still fit; only the least intense isotope
  labels can differ when there is no room for them.
* Shading of the time ranges on the chromatograms: on by default (as in the
  window); untick it in the report window to leave it out.
* Copied and saved images of a zero charge mass spectrum: the labels are
  placed for the size of the image (they ran into each other in a journal
  sized image), the description above the plot is set smaller or on two
  lines (it pushed the plot aside), and the x axis has as many numbers as
  fit (they overlapped).
* Report headline: a label too long for its box is written on two lines (it
  ran into the next box).
* Updates: a library in use while an update is copied (msengine.dll) can be
  delivered as msengine.dll.new; it replaces the old one at the next start.

3.1 (03.10.2026)
* UniDec runs inside MSpektra: its engine is built into the program's
  library (msengine.dll, compiled from UniDec's own C source) instead of
  being started as the separate program unidec.exe. Same results (masses,
  intensities, scores and R squared agree with unidec.exe to the last
  printed digit, the spectra to single precision rounding), 1.4 to 4.5 times
  faster on the test spectra, about 70 to 90 MB less memory, no temporary
  files, and Cancel stops a run within a fraction of a second. The notes of
  a result end with "UniDec engine run in the library (n threads)". The
  engine is UniDec's own (Michael T. Marty, Marty et al., Anal. Chem. 2015);
  unidec.exe stays in the package for the classic UniDec windows of the
  start screen.
* Settings where unidec.exe crashed now give a plain message, and UniDec
  results no longer depend on the number of processor cores (two places in
  the original engine where its threads raced, and a few reads past the end
  of an array, are fixed).
* The remaining calculations run in the library as well: the calibrant
  search and the calibration models (fits, cross validation, applying them),
  the apex of profile peaks, isotope patterns, the formula finder and the
  isotope match, LC peak integration and smoothing, PDA chromatograms and
  UV spectra, the peak labels of spectra, and the steps around each
  deconvolution (data preparation, minimum intensity, baseline, isotope
  labels, most abundant isotope). Checked against the previous Python code
  on the test files: identical, apart from rounding in the last digits.
  Faster: the isotope pattern of a 36 kDa protein from 3.5 s to 0.12 s, the
  calibration model table from 65 ms to 2 ms, the formula finder from 39 ms
  to 4 ms, the calibrant search from 0.47 s to 0.03 s; the first peak
  integration of a session no longer loads scipy (about 3 s).
* Smoother: the step that thins a large spectrum for drawing takes under
  0.01 s instead of about 0.2 s for 600 000 points (the same points are
  drawn), so zooming, panning, clicking a scan and averaging redraw sooner;
  the start screen appears sooner (NumPy loads only when it is needed); the
  text files written after each deconvolution take about a third of the
  time (same files).
* Label all peaks in view: a selected, pinned or added mass keeps one label
  with its own mass (its peak was labelled a second time); stacked labels
  stay inside the plot; IsoDec results are labelled too.
* Reports of isotope resolved results: a table of the species (average mass,
  most abundant isotope, isotope peaks, intensity, share) and the isotope
  peaks of the main species (of every species when expected masses are
  entered), instead of one M+0, M+1, ... list across all species; the labels
  of the mass figure and the charge states of the spectrum figure no longer
  overlap; masses added by hand are marked "(added)"; the most abundant
  isotope is given with the next one when the two are equally tall, as in
  the window.
* Files with positive and negative scans: calibrating the second polarity
  removed the calibration of the first. Each polarity now keeps its own, the
  Calibration chip shows both, and Remove in the calibration window removes
  only that polarity's calibration.
* The Deconvolute and calibration windows of a file that was already open
  offered old settings (and the calibration panel could write them back over
  newer ones); they now offer the last settings used.
* LCMS Full window button: negative ion spectra opened as positive ions.
* Several Bruker .d folders or mzML files dropped onto MSpektra.exe now
  open together (only the first opened).
* Masses added by hand to an isotope resolved result are labelled with their
  mass only (no "average" line), and the table shows "(added)" in full.
* When every result is hidden, the report window says to show one again.
* MSpektra.exe is the only program file in the folder: the old launcher
  UniDec.exe of the UniDec named versions (it started the same program) is
  moved to backups\old_launcher at the first start.

3.0 (03.10.2026)
* Isotope resolved species: the most abundant isotope is now chosen from the
  isotope peaks in the m/z data (their areas at every charge state, summed),
  no longer from the heights after deconvolution, which can tip two nearly
  equal isotopes by a percent. The table and the tooltip give the second
  isotope with one decimal (e.g. "12228.2451 (12229.2485: 99.7 %)"); the
  tooltip says that the ratio comes from the m/z data. On Protein POS.d,
  3.17 to 3.31 min, the two isotopes are 100 and 99.7 % in the data itself
  (per charge state the ratio scatters from 95 to 109 %), so both are named;
  over 3.21 to 3.29 min 12229.250 is the taller one by 2 %.

2.9 (03.10.2026)
* Isotope resolved results: when the isotope next to the most abundant one
  is at least 95 % as tall, both are named: the label of the species reads
  e.g. "12228.2451 / 12229.2485" with "average 12228.99, isotopes within 1 %"
  above it, and the table gives "12228.2451 (12229.2485: 99 %)". For a
  protein of about 12 kDa the isotopes +7 and +8 are nearly equal (averagine:
  100 and 95 %), so which one is the tallest depends on the noise and on the
  scans averaged (Protein POS.d: 3.17 to 3.31 min gives 12228.245 at 100 %
  and 12229.249 at 99 %; 3.21 to 3.29 min gives 12229.250 at 100 % and
  12228.246 at 97 %). The masses themselves were right; one of two equal
  peaks was named.
* Peak labels in a zoomed mass spectrum: they were placed for the size of a
  copied (journal sized) image, so a zoom of 20 Da on screen got only the
  species label. They are now placed for the plot on screen: the isotope
  peaks in view are labelled (up to the number of labels set in the
  dialog), and a wide view labels the species instead of stacking isotope
  labels in one column.
* Right click on a mass spectrum > Label all peaks in view: every peak of the
  current view gets its mass (kept while zooming; small ghost bumps half way
  between isotope peaks are left out); "Automatic peak labels" goes back.

2.8 (02.10.2026)
Deconvolution models (maximum entropy and UniDec, LCMS and HRMS)
* Maximum entropy, rebuilt on the published method (Skilling and Bryan;
  Ferrige et al. for electrospray): the fit now uses a noise model measured
  on the spectrum itself (a noise floor plus counting noise that grows with
  the intensity; no settings needed), so weak signal is fitted as carefully
  as strong peaks and the strong peaks' shape no longer creates ghost masses
  (harmonics at M x z'/z: B Lac 8 extra peaks in 2.7, 1 now). The charge
  state envelope may differ between species: after a first pass the program
  places envelope nodes at the species it found, so mixtures of proteins with
  different charge state distributions come out right; a species never gets
  one charge state alone. The weight of the entropy is searched for the
  classic criterion (chi squared per point 1); on real data the peak shape
  keeps the best fit above it and the notes then say "target 1 not reachable
  with this peak shape". The notes list the model used (noise floor, counts,
  envelope nodes). It converges in about half the rounds of 2.7.
* Maximum entropy, wide mass ranges at fine steps: the zero charge axis was
  capped at 400000 points, so 1000 to 50000 Da every 0.01 Da had a point
  every 0.12 Da: isotope peaks were drawn from two points and the marked
  apex sat off the drawn maximum. The axis now keeps the mass step (up to 6
  million points, thinned for display away from the peaks). On the 12 kDa
  protein (3.21 to 3.29 min, 1000 to 50000 Da, 0.01 Da, resolved): 2.7 found
  one mass (most abundant isotope 12229.278, R squared 0.87, 6 min); 2.8
  finds 12229.249 and the oxidised form 12245.243 at 22 % (R squared 0.99,
  3.5 min on 2 cores).
* UniDec: the crash "the deconvolution engine stopped with code 3221226505"
  is fixed. It came from the UniDec engine itself: when no peak of the m/z
  range has the same mass at the next charge state, the engine removes every
  point and stops (Windows code 0xC0000409). This happened without data
  reduction when a minimum intensity left only isolated peaks (for example
  singly charged ions), and always with Charge smoothing 0. Such data are now
  put on bins of their own spacing (the notes say so); settings the engine
  cannot run are refused before it starts, with the reason. Engine failures
  are shown in plain words, and a failed run is no longer repeated.
* UniDec result quality: every result is rated good, fair or poor from
  UniDec's own scores (the DScore of every mass and the UniScore, Kostelic
  and Marty 2022), with one or two sentences on what is wrong and what to
  try. A spectrum without a charge state envelope in the m/z range (for
  example only singly charged peaks) gives hundreds of "masses" that are one
  peak read at many charge states; such a result is now marked "Poor result"
  in red on its row, above the mass table and in the list of files, and its
  notes start with the explanation. The mass table has a Score column.
* UniDec, isotope resolved (mass steps below 0.1 Da): the peak width is now
  right over the whole m/z range (12 kDa protein: R squared from 0.73 to
  0.95), and weak species no longer take in their neighbours (average masses
  within 0.1 Da).
* UniDec memory check: from what the engine really allocates (measured);
  settings are refused only when they would not fit, with the need in the
  message. Mass steps finer than the engine's 32 bit mass axis can hold are
  refused (for example 0.001 Da at 14000 Da). The Isotope mode setting has
  no effect in UniDec 8.2.1 (its engine ignores it); the notes no longer
  claim monoisotopic masses.
Deconvolution results in the windows
* The charge states plot is off by default; right click a result > Show
  charge states (Hide charge states removes it). The choice is remembered.
  The mass spectrum then takes the whole width.
* The list of open files shows the deconvolution results of each file as
  branches, newest first: click one to show it in the table and scroll to
  its row, the eye hides or shows it (kept), x closes it. Keys after a click
  in the list: Up/Down, Enter, Space (hide/show), Delete (close), Left/Right
  (fold). A red triangle marks a poor result.
* The shaded area under the mass spectrum was closed by a diagonal from the
  bottom right to the first visible point (a sloped band when the view began
  on a peak); fixed.
* Isotope resolved results: the dot of the selected species sits on the top
  of its most abundant isotope peak; its bold label is that isotope's mass,
  with the average mass in small type above it. The table has "Average mass"
  and "Most abundant isotope". Maximum entropy results without resolved
  isotopes are no longer listed with those columns.
* Presets now also set the mass range; the maximum entropy fit is drawn on
  the whole spectrum; labels stay inside a tile made shorter; zoomed mass
  axes show whole numbers instead of an offset.
HRMS Analysis
* The calibrant segment is found again when a file opens (it was not with
  the C++ core of 2.4 to 2.7) and "Calibrate automatically" works.
* The calibration window uses every scan of the calibrant range (it dropped
  the first and last); in files with both polarities a calibration applies to
  its own polarity only (the chip says "(+)" or "(-)").
* Spectrum titles show the right scan time in files with polarity switching
  (also LCMS Analysis).
* Reports show sample name, instrument, operator, method and acquisition
  date of Bruker files and mention the DataAnalysis recalibration; the model
  table matches the calibration applied. The recalibration note shows the
  latest DataAnalysis calibration.
* Labels, the formula check and the report give the same m/z for a peak.
  The formula finder with [M]+ and [M]- also finds even electron ions;
  formulas such as CuSO4.5H2O are read correctly; the formula check takes the
  peak nearest to the calculated m/z. HPC calibration uses at most n - 8
  terms on n calibrants.
* Spectrum exports (text, JCAMP-DX, clipboard) write m/z with 6 decimals.
* mzML: integer arrays and times given in seconds by accession are read;
  MS-Numpress files are refused with a message on how to convert them; a
  scan that can no longer be read is reported.
LCMS Analysis
* Files with positive and negative scans are split by the scan event stored
  with every scan, so a single damaged scan no longer mixes the polarities.
* Saturated scans are marked again with the C++ core (red ticks and the
  "saturated signal left out" note).
* Bin width: 0.001 to 10 m/z; anything else is replaced by the last valid
  value with a note. Smoothing "3 points" now smooths. A wavelength change in
  the PDA view keeps an averaged UV spectrum. A new default window for mass
  chromatograms applies to new ones only. Ranges and backgrounds outside the
  run, exports that cannot be written (file open in another program) and
  unreadable files are explained in plain words. Several .lcd files dropped
  on the program open together.
Robustness
* The C++ core was run under memory and thread sanitizers with 227 test
  cases (damaged mzML, .lcd and .d files, wrong arguments, threads,
  cancellation): 75 problems found and fixed (a memory error on every .lcd
  open, crashes when a file was closed while it was being averaged, two
  threads using one file, endless calculations in the formula tools, a wrong
  isotope pattern above 215 kDa on Windows). Damaged or truncated files are
  reported instead of stopping the program.
* Settings are saved safely when several windows of the program save at the
  same time; the deconvolution log records each job and, after a crash,
  where it happened.

2.7 (02.10.2026)
* The fixed size limits of the deconvolutions are gone: no more "at most
  400 000 mass points" for maximum entropy and no more "1 000 000 mass
  points / 60 million cells" for UniDec. Both now ask the library how much
  memory this computer has (half of the physical memory is allowed for one
  deconvolution) and refuse only settings that would not fit, with the
  estimated need in the message. The maximum entropy model takes less than
  half the memory it did (10 bytes per entry instead of 24; the envelope
  weights are applied on the fly) and builds without a temporary copy.
  Tested on the 12 kDa protein, resolved isotopes at 0.01 Da (2 cores):
  10 000 to 16 000 Da (600 000 points) 82 s and 0.5 GB; 5 000 to 25 000 Da
  (2 million points) 4.3 min and 1.3 GB; 1 000 to 50 000 Da (4.9 million
  points) 5.6 min and 1.5 GB. UniDec at 0.01 Da over 1 000 to
  50 000 Da (4.9 million points) 12 s, at 0.001 Da over 11 000 to 14 000 Da
  (3 million points) 28 s; the engine crash that led to the 1 000 000 point
  guard in 2.3 came from unbinned data with gaps (fixed there), not from
  the number of mass points. The Python code (used only when the library
  is turned off) keeps the old limits.
* The notes of a maximum entropy result from the library no longer repeat
  the m/z range and minimum intensity line.
* Wide mass ranges at a fine step give a lot of mass points to draw and to
  save (a 3 million point result is a 90 MB _mass.txt): the window stays
  usable, but prefer a coarser step unless resolved isotopes are needed.

2.6 (02.10.2026)
* Shimadzu files: the bin width of the LCMS window (any value, not only
  the default 0.05 Da) and the background subtraction now run in the
  library as well; the Python reader is no longer loaded for them. All
  combinations of bin width, background ranges and fill checked against
  the Python code on the two test files: identical (the 2.5 background
  subtraction had differed below 0.1 % of the base peak).

2.5 (02.10.2026)
* Second step of the conversion to C++ (msengine.dll): LCMS Analysis now
  reads Shimadzu .lcd files (MS scans, scan events, PDA data) through the
  library; the UniDec deconvolution (preparation of the data, the UniDec
  engine, peak picking, isotope grouping) and IsoDec (through UniDec's
  isodeclib.dll) run in C++ as well, for all three windows. Checked against
  the Python code on the test files: Shimadzu chromatograms, averaged and
  single spectra, PDA data and sample information are identical; UniDec
  masses, peaks, fit and R squared are identical (the UniScore differs in
  the fourth digit); IsoDec results are identical. Times (2 cores): UniDec
  of the 12 kDa protein 43 s -> 7 s including the UniDec engine, IsoDec
  13 s -> 1 s, opening a 100 MB .lcd 1.4 s -> 0.3 s. Every result from
  the library says "C++ core" in its notes.
* The isotope labels of resolved UniDec results no longer import scipy's
  signal module (4 s on the first result in a fresh deconvolution process).
* Also in the library (for the C++ interface of the next step): the LC
  peak integration, the calibration models, the formula tools (mass,
  isotope pattern, formula finder) and the Shimadzu sample information.

2.4 (02.10.2026)
* First step of the conversion to C++: the computational core of HRMS
  Analysis is now a C++ library, _portable\msengine\msengine.dll (source in
  the cpp folder of the project; built with GCC/MinGW, OpenMP, one self
  contained DLL): reading Bruker .d folders (analysis.baf, analysis.tsf,
  analysis.tdf through Bruker's own libraries) and mzML files, TIC, base
  peak and mass chromatograms, averaging with background subtraction,
  centroiding, peak apex, and the maximum entropy deconvolution (LCMS and
  HRMS). The windows, plots and reports are still Python (next steps).
  Results were checked against the Python code on the test files: spectra
  and chromatograms are identical to the last digit; the maximum entropy
  uses its own optimizer, masses agree within 0.03 Da and heights within
  10 %. Times (same computer, 2 cores): mass chromatogram of a 100 MB
  mzML 1.7 s -> 0.02 s, maximum entropy of the 12 kDa protein (resolved
  isotopes) 30 s -> 7 s, wide range 1000 to 50000 Da 7 s -> 1 s, low
  resolution LC-MS 14 s -> 4 s. Each result says "C++ core" in its notes.
  When the library cannot be loaded, the Python code is used (a line in
  the log says so); "engine": "python" in config\portable_settings.json
  turns the library off.
* OpenMP threads of the library do not spin wait (OMP_WAIT_POLICY passive),
  so the interface stays responsive on computers with few cores.

2.3 (02.10.2026)
* Start screen: the Deconvolute window (and the classic UniDec tools) now
  start as their own process, so the start screen and LCMS/HRMS Analysis no
  longer load the UniDec interface libraries in the background. A click
  on a tile no longer waits for that loading (it waited up to the whole
  background load, 10 to 16 s on a cold start). The background loading of
  the libraries LCMS/HRMS need is shorter and no longer blocks a click; a
  tile ignores further clicks while its window is being built (a second
  click during the loading crashed the app). Setting
  "deconvolute_own_process": false in config\portable_settings.json
  restores the old in-process Deconvolute window. _portable ships its
  compiled Python files; MSpektra.bat no longer runs the PowerShell
  Unblock-File pass over all files on the first start (the launcher does
  the needed part itself, once).
  Note on "multithreading": importing libraries in several threads does
  not make a start faster (Python runs one thread at a time for this),
  it only slows the window being built. The deconvolution methods already
  run in a separate process (all cores for UniDec and, from 2.3, maximum
  entropy), and the Deconvolute window is now a process of its own.
* LCMS and HRMS: every deconvolution adds a new result row (zero charge
  mass spectrum and charge states) above the older ones; the newest is at
  the top, right below the spectrum. The mass table and the fit on the
  spectrum belong to the active result (the newest, or the one whose plot
  was clicked last; its row says "shown in the table"). Right click a
  result > Hide this result / Hide all results. At most 8 results are
  kept (the oldest is dropped). Each row's subtitle names the method,
  time range and mass range.
* HRMS Analysis: the deconvolution result rows open at full size while the
  chromatogram and the spectrum stay at half size.
* Deconvolute dialog: switching the method (or a preset or the polarity)
  keeps the position and the size of the dialog.
* Maximum entropy: a mass step too coarse for resolved isotope peaks (e.g.
  0.5 Da at 12 kDa, where an isotope peak is 0.2 Da wide) made the fit
  impossible, and with a wide mass range the result was dozens of masses
  of similar height at fractions of the real mass (M x z'/z). Such a
  step now fits isotope envelopes (average masses) instead, with a note
  in the result; use 0.01 to 0.05 Da and a narrow mass range for resolved
  isotopes. The calculation uses all processor cores (one is left free on
  computers with more than 4) and a sparse model: 1000 to 50000 Da on the
  12 kDa protein 192 s -> 8 s, the resolved 11000 to 13500 Da case 85 s
  -> 30 s (same results). The mass axis of wide ranges is thinned away
  from the peaks (the _mass.txt file was 12 MB).
* UniDec: more than 1000000 mass points (mass range / step) are refused
  with a message: the engine could stop with a crash (code 0xC0000409)
  above that, e.g. 11000 to 14000 Da at 0.002 Da.
* HRMS Analysis: averaging a long time range of a Bruker .baf file whose
  profile spectra have different m/z axes (zero trimmed) collected all
  points of all scans first: 1000 scans took gigabytes of memory, the
  computer swapped and the window froze for minutes. The scans are now
  merged into bins one by one (memory of one spectrum). An mzML file with
  profile spectra kept every scan in memory (a 10 min run at 1 Hz: 5 GB);
  now at most 400 MB are kept and the other scans are read again from an
  indexed mzML when needed. Long averages (HRMS) run with a progress
  window and Cancel, and the window stays responsive.
* Fixes from a code review: formula check and peak labels on centroid
  only data (mzML with peak picking) took a Gaussian fit over neighbouring
  sticks (hundreds of ppm off); the report crashed on an IsoDec result
  with no masses and on a chromatogram whose peaks have zero area; IsoDec
  results with masses spread over a wide range were drawn on millions of
  points (slow, 150 MB files); a failed saving of a result left the panel
  busy; settings edited by hand with a wrong type stopped the Deconvolute
  button; the number of labels can be 0; the formula finder window is
  closed with its file; the temporary copy of a write protected .d folder
  is named after the file (two files of the same name shared one copy);
  the worker log is started afresh above 2 MB; cancelled worker processes
  left pipe handles open; the settings file is written atomically.
* Everything else is unchanged from 2.2.

2.2 (01.10.2026)
* Maximum entropy with a wide mass range gave masses that do not exist
  ("ghosts" at simple fractions of the real mass, e.g. 10/9, 2/3, 5/8 of
  it), often taller than the real one. Cause: the starting charge
  envelope came from a least squares fit on a flat mass spectrum, which
  puts all weight on a few charges; the masses then settled on a single
  charge, each charge state of the protein becoming a mass of its own.
  The envelope now starts flat over every charge with signal, and for 15
  charges or more it is kept smooth. The real mass is now the tallest
  for mass ranges such as 1000 to 50000 Da (tested on an intact protein
  and on LC-MS spectra).
* UniDec: the masses are the apex between the points of the mass grid,
  not the grid point itself. With a 0.1 Da mass step the grid alone put
  a 12 kDa mass up to 4 ppm off (12229.20 instead of 12229.25 Da in the
  test file). With resolved isotopes (mass step up to 0.25 Da) the table
  now has one row per species (average mass, most abundant isotope,
  number of isotope peaks), as for maximum entropy, instead of one row
  per isotope peak. A peak window below two mass steps is raised to two
  steps (it made every point a peak).
* Isotope groups of two species next to each other (e.g. +16 Da) are no
  longer merged into one.
* HRMS Analysis: m/z labels of the spectrum (and the base peak read out,
  and the report) are at the apex of the profile peak, not at its
  highest data point (on a TOF spectrum the points are 3 to 5 ppm apart,
  so labels were up to 2 ppm off). They now agree with Bruker's peak list
  within about 0.3 ppm. Clicked peaks and the calibration already used
  the apex.
* HRMS Analysis: a recalibration saved by DataAnalysis in the .d folder
  (calibration.sqlite) is used, and the Calibration chip and the file
  summary now say so ("file recal. (DataAnalysis, date)"). A calibration
  made here comes on top of it.
* Everything else is unchanged from 2.1.

2.1 (01.10.2026)
* UniDec failed with "No such file or directory" for data files with long
  names in deep folders (e.g. OneDrive): UniDec writes
  <name>_unidecfiles\<name>_..., which went past the Windows limit of 260
  characters for a path. When that would happen, UniDec now works under
  a shortened name (start of the name plus a short code); the results
  (_mass.txt, _peaks.csv, input) keep the full name. The same applies to
  the spectrum handed to the Deconvolute window (right click > open in the
  Deconvolute window). Maximum entropy and IsoDec were not affected.
* Deconvoluted mass spectrum (LCMS and HRMS): right click a peak > Label
  ... (keep the label) keeps its mass label at every zoom and in the
  report; right click where no mass was found > Add the mass here adds
  the highest point there to the masses (label, table marked "added",
  _peaks.csv, report); right click it again > Remove the added mass.
* HRMS Analysis: tiles open at half size by default (right click > Tile
  heights still offers "first 2 tiles fill the window").
* Everything else is unchanged from 2.0.

2.0 (01.10.2026)
* HRMS Analysis reads timsTOF data: .d folders with analysis.tdf and
  analysis.tdf_bin. Each MS1 frame is one spectrum with the ion mobility
  dimension summed (the SDK's profile of the frame); MS/MS frames (PASEF,
  DIA) are counted but not shown. TIC, base peak (the frame maxima
  recorded by the instrument), mass chromatograms (extracted by the SDK
  when asked for, so opening is immediate), averaging, background,
  calibration, deconvolution and reports work as for the other formats.
  Same TDF SDK library as in 1.9. Ion mobility itself (1/K0 filtering,
  mobilograms, CCS) is not shown. Everything else is unchanged from 1.9.

1.9 (01.10.2026)
* HRMS Analysis reads the newer Bruker format: .d folders with analysis.tsf
  and analysis.tsf_bin instead of analysis.baf, as written by otofControl
  6 and later (e.g. maXis II). Profile and line spectra, TIC, base peak
  and mass chromatograms, averaging, background, calibration, formula
  check, deconvolution and reports work as for analysis.baf. Read with
  Bruker's TDF SDK (timsdata.dll, in _portable\timsdata with its licence
  files); the .d folder is only read, nothing is written into it.
  timsTOF data (analysis.tdf) were not yet supported (see 2.0).
  Everything else is unchanged from 1.8.

1.8 (30.09.2026)
* LCMS and HRMS: the deconvoluted mass plot uses the requested mass range
  for its initial view and zoom-out limit, even when the engine returns a
  narrower data range. Full view restores the complete requested range.
  Applies to UniDec, Maximum entropy and IsoDec. Empty results keep the
  requested axis; single-mass results have a usable plotting range.
* Charge plots retain the requested charge range. Zoom to masses found
  handles peaks in any order and remains reversible with Full view.
* Includes the deconvolution precision (code 104), validation, worker,
  export and HRMS averaging fixes made earlier on 30 September.

1.6 (30.09.2026)
* Reports: the mass spectra (and the UV spectrum) are the ones shown in
  LCMS Analysis (same time range, background and bin width); the caption
  says which. Before, they were averaged again over the main peak without
  the background subtraction, so they could differ from the window.
* Deconvolution: the ranges typed in are kept. "Suggest them every time
  this window opens" is now off by default (it was on in HRMS Analysis and
  replaced the mass range each time). After a deconvolution the zero
  charge spectrum shows the whole mass range deconvoluted; right click >
  Zoom to the masses found narrows it.
* Calibration window: the selected calibrant (clicked in a plot or the
  table) is ringed in orange in the spectrum and the mass error plot, with
  a dashed line at its m/z.
1.5 (30.09.2026)
* Smooth zooming, panning and dragging in the zero charge mass spectrum
  (and the other spectra): dense traces are drawn from their highest and
  lowest point per pixel column (identical on screen, several times
  faster; saved images use their own resolution), the drag box is drawn
  over a stored image of the plot, and after a zoom only the mass labels
  are drawn again.
* Reports: "Figures in the report" in the report window, a tick box for
  every figure (chromatogram, UV spectrum, mass spectrum and total ion
  chromatogram of each polarity, isotope patterns, deconvolution,
  calibration). The choice is kept. The shaded time range on the total
  ion chromatograms is now off unless ticked; "(slider)" is no longer
  printed after the wavelength.
* LC-MS purity report: every graph at full width and the same height,
  larger than before: chromatogram, total ion chromatograms, UV and mass
  spectra of the main peak, then extracted ion chromatograms (those set in
  LCMS Analysis with Mass chrom., or else the base peak of the main peak).
* HRMS Analysis: the chromatogram and the spectrum fill the window; further
  tiles follow below (scroll). Right click a tile > Tile heights for the
  other sizes.
* Calibration window: opens almost full screen; the calibrant spectrum and
  the mass error plot fill it (the calibrant table below). Right click a
  point: "Remove this point from the fit"; or click a point and press
  Delete.
1.4 (30.09.2026)
* Faster start: the start screen needs only the window library and appears
  in about 1 to 2 s; the libraries of the windows load in the background
  meanwhile. LCMS and HRMS Analysis no longer load UniDec's interface at all;
  the Deconvolute window loads it when it is opened. The deconvolution
  process starts in the background after an analysis window opens, so the
  first deconvolution does not wait for its libraries. The first-start
  preparation checks only program files (DLL, EXE), not all 11 000 files.
* Deconvolution window: m/z range and minimum intensity in "Spectrum and
  data", with the time range, base peak and the peaks and points used.
  The minimum intensity is in the units of the spectrum's y axis (counts)
  or in % of the base peak; peaks whose top is below it are removed and
  the peaks above it are kept whole; it applies to every method. It
  replaces UniDec's own intensity threshold, which worked on normalised
  data (0 to 1) and so removed everything with values such as 40.
  Fixed: UniDec stopped with "int() argument ... NoneType" when nothing
  was left above that threshold. UniDec refuses a single charge state
  (its engine stops with one; use 1 to 2, or maximum entropy).
* The deconvolution progress window is centred on the window.
* LCMS Analysis, HRMS Analysis and the Deconvolute window open maximised.
* Plots: black axes, ticks and labels, blue traces, thinner lines (also in
  exported images, reports and the Deconvolute window).
* Sharper start screen and custom controls: their text is drawn by
  Windows' own text renderer (ClearType) in Segoe UI below 13 pt; icons
  sit on whole pixels.
* Timings in the log: "Ready for the first window", "start screen shown",
  the background loading, "window built", "read <file> in".
1.3 (29.09.2026)
* Negative ions throughout: the deconvolution window follows the polarity
  of the spectrum and offers the negative charge carriers (loss of H+,
  Cl-, formate, acetate, custom); UniDec, maximum entropy and IsoDec use
  them; charge states are labelled z-. Fixed: UniDec stopped with an error
  on negative spectra (its automatic peak width came out negative).
* Check a formula starts with an ion of the spectrum's polarity ([M-H]-
  for negative ions) and says so when ion and spectrum do not match.
* Negative spectra sent to the Deconvolute window open in negative mode;
  a .lcd file opened in the Deconvolute window opens in LCMS Analysis.
1.2 (29.09.2026)
* Reports: HRMS compound report, LC-MS purity report, deconvolution report
  and Supporting Information page, as PDF or Word; the kind is chosen each
  time (Report... in the top bar, Ctrl+R, or right click any plot).
* Several files open in one window: file panel on the left (open files,
  switch with a click or Ctrl+Tab, close with x or Ctrl+W without closing
  the window) and, as in LabSolutions, the data files of the folder with
  date, size and sample name (double click to open).
* Sample information of .lcd files (sample name and ID, acquired by and
  when, vial, injection volume, method) in the file panel and the reports.
* Base peak (m/z and intensity) in the title of every spectrum.
1.1 (29.09.2026)
* Deconvolution window with every setting of UniDec, maximum entropy and
  IsoDec; runs in a separate process with Cancel; settings that cannot
  finish are refused. Maximum entropy reworked (matches DataAnalysis
  MaxEnt on G4COL within 1.3 ppm); absolute intensities everywhere.
* Calibration: remove calibrants from the fit; scrolling tiles.
* Tiles: half size by default, each tile resizable on its own, Close tile.
* Images: copied and saved at 8.5 x 6 cm, 9 pt Arial, 300 dpi, margins
  fitted to the labels; copying uses the Windows clipboard directly;
  Ctrl+C / Ctrl+S on a plot; a note confirms each copy or save.
* Plots: ticks outside the frame, fewer y ticks (they follow the tile
  height), room above the tallest peak for its label.
* LCMS: MS and PDA times linked; integrated peaks kept (areas measured
  again) when the wavelength or trace changes; marker lines left out of
  images.
1.0 (29.09.2026)  first release of MSpektra (LCMS Analysis, HRMS
  Analysis, Deconvolute window on UniDec 8.2.1).


QUICK START
-----------
1. BEFORE extracting: right-click the .zip > Properties. If there is an
   "Unblock" checkbox at the bottom, tick it and click OK.
2. Extract to a local folder with a short, plain path, for example
      C:\Users\<you>\MS_Analysis      or      D:\Tools\MS_Analysis
   Best is any folder you can copy into WITHOUT an administrator prompt
   (your user folder, Documents, a USB stick). Avoid OneDrive/SharePoint-
   synced folders, network drives and paths with special characters.
   If you do put it in a protected place (C:\Program Files), it still
   works: settings, recent files, logs and caches then go to
   %LOCALAPPDATA%\UniDecPortable instead of this folder.
3. Double-click  MSpektra.exe . The start screen has three tiles:
   LCMS Analysis, HRMS Analysis and Deconvolute. (The classic tools of the
   UniDec package, MetaUniDec, UniChrom, IsoDec, Data Collector, ..., are
   hidden: right click the empty background of the start screen.)
   MSpektra.bat does the same and stays as a fallback in case security
   software blocks MSpektra.exe.
4. Taskbar / Start menu: double-click  Create_desktop_shortcut.bat  once. It
   makes "MSpektra" shortcuts on the desktop and in the Start menu and
   removes older "UniDec" shortcuts. To pin: start MSpektra, right-click
   its icon in the taskbar, choose "Pin to taskbar".
   Run Create_desktop_shortcut.bat again if you move the folder.
   Updating from the UniDec-named version: the old launcher UniDec.exe is
   moved to backups\old_launcher at the first start (it is not used);
   UniDec.bat and UniDec_console.bat can be deleted after this step (an old
   pinned UniDec icon should be unpinned and MSpektra pinned instead).

Opening a data file directly: drag it onto MSpektra.exe, or use
"Open with > MSpektra.exe": LC-MS and HPLC files (.lcd, Agilent .D, Waters
.raw, Thermo .raw, ANDI .cdf) open in LCMS Analysis; a Bruker .d folder (or
analysis.baf inside it), an Agilent MassHunter .d, a Thermo Orbitrap .raw
or an .mzML file in HRMS Analysis; anything else (.jdx, .txt, ...) in the
Deconvolute window.

MSpektra.exe is a small launcher that starts the bundled Python; it is
not digitally signed, so a copy downloaded from the internet may show a
one-time Windows "unknown publisher" prompt.


LCMS ANALYSIS: SHIMADZU .LCD FILES
----------------------------------
LCMS Analysis reads LabSolutions LC-MS data files (.lcd, e.g. LCMS-2020 with
PDA) directly: no export from LabSolutions is needed. Start it from the
start screen, or drag a .lcd file onto MSpektra.exe.

Right click first (LCMS Analysis and HRMS Analysis)
* Everything is done from the plots: right click a chromatogram or a
  spectrum and choose; operations that need settings open a small window
  (average a time range and background, mass chromatogram with its window,
  integration, deconvolution, calibration, formula check, PDA settings).
* The side panel is secondary and hidden by default: "Show panel" / "Hide
  panel" in the top bar (or F9) shows or hides it, and the choice is
  remembered. Nothing needs it.
* Mouse: Ctrl + drag zooms a box on any plot (with any tool), double click
  gives the full view. The mouse wheel scrolls the column of tiles
  (Ctrl + wheel zooms the axis).
* Tiles open at half the height of the journal proportions (a copied image
  is 8.5 x 6 cm), so several plots fit in the window; the column scrolls
  when it is longer than the window (scroll bar on the right, or the mouse
  wheel). Drag the gap below a tile to make only that tile taller or
  shorter; the others keep their height. Right click > Tile heights: "Half
  size" (default), "Journal proportions (as copied images)" or "Fit every
  tile in the window" (double click a gap does the last). Tables keep their
  own height. The same applies in the calibration window. Copied and saved
  images do not depend on the tile size.
* Open files (left panel): every file opened (Open, Ctrl+O, drag and drop,
  or a double click in the folder list) gets its own views and keeps its
  ranges, peaks, calibration and results. Click a file to show it
  (Ctrl+Tab for the next one); x on the file, Ctrl+W or right click >
  Close closes it without closing the window (also "Close the other
  files", "Close every file"). A file that is open already is just shown.
  Only one deconvolution runs at a time (for all files).
* Folder (below the open files): the data files of the folder of the file
  shown (or of the last folder), newest first, with date, size and sample
  name, as the data browser of LabSolutions. Double click opens a file;
  the tooltip shows the sample information; the folder icon shows another
  folder; the list follows new runs by itself. "<" narrows the panel.
* Right click any tile > "Close this tile" hides it (e.g. the charge states
  or the UV spectrum); "Closed tiles" in the same menu shows it again.
* Intensities are absolute everywhere (counts as recorded, mAU for the PDA),
  including the deconvoluted spectrum. "Scale each trace to 100 %" (right
  click a chromatogram) is still there when you want relative traces.

Mass spectrometry view (everything for the mass spectra in one view)
* One tile per scan event: positive and negative scans (polarity switching)
  are split automatically and shown side by side, each with its chromatogram
  (TIC or base peak) and its mass spectrum. The time axes are linked.
* Drag across a chromatogram to average that time range (both polarities);
  Shift + drag sets a background range that is subtracted. A click shows the
  spectra of the single scan at that time (marked by a dashed line; the
  range highlight is cleared).
* Mass chromatograms: right click a chromatogram > "Add a mass
  chromatogram" (m/z, window in m/z or ppm, scan event), right click a peak
  in a spectrum > "Mass chromatogram of m/z ..." (default window) or "...
  with options" (the window chosen there becomes the default), or type m/z
  lists in the side panel, each with its own window if wanted: 803.2,
  1204.8+-0.3, 889.09+-10ppm. Right click a mass
  chromatogram > "Edit this mass chromatogram" changes m/z and window, or
  removes it. They appear as their own tiles below the TIC ("Stacked
  tiles"), or on top of the TIC ("Overlaid on the TIC").
* Spectra are shown as a profile on the full m/z grid (the line returns to
  zero between peaks), or as sticks (Spectrum display).
* Right click a spectrum > "Deconvolute this spectrum" (or Deconvolute in
  the top bar) opens the deconvolution window with every setting; the
  results appear in the same view (see DECONVOLUTION). LCMS Analysis also
  offers "Open this spectrum in the Deconvolute window" (the UniDec
  program); HRMS Analysis does not send spectra there, everything is in its
  own window. Export saves .txt or .jdx.
* Right click a chromatogram also offers: Trace (TIC, BPC, mass
  chromatograms only), how mass chromatograms are shown, Smoothing, Scale
  each trace to 100 %, Link the time axes, Clear the background range.
* Everything saved for a data file (spectra, chromatograms, peak tables,
  images, deconvolution results) goes into one folder next to it:
  <file name>_analysis.
* Saturated detector readings are left out, as in LabSolutions; the scans
  concerned are marked by small red ticks at the top of the chromatogram.

PDA tab
* Wavelength map (time x wavelength) at the top with two sliders: the round
  handle on the top edge sets the time, the handle on the right edge sets
  the wavelength. Drag the handles or the crosshair, or use the arrow keys.
  The chromatogram at that wavelength and the UV spectrum at that time
  follow immediately.
* Chromatograms use the bandwidth set (as "280nm,4nm" in LabSolutions);
  more wavelengths, max plot, smoothing and an overlay of the MS trace with
  an adjustable detector delay are available. UV spectra can be averaged
  over a dragged range, with background subtraction.
* Integrated peaks stay when the wavelength slider moves (or the bandwidth,
  smoothing or trace changes): each peak keeps its start and end time and
  its area and height are measured again on the new chromatogram. The same
  holds for MS traces (e.g. a mass chromatogram whose window is edited).
* MS and PDA times are linked: a click or double click on the PDA
  chromatogram or map sets the same time in the MS view (the MS spectrum of
  that peak is ready when you switch back), and a click or a dragged range
  on the MS chromatogram moves the PDA time slider. The detector delay set
  in the PDA settings is taken into account. Right click a chromatogram >
  "Link MS and PDA times" switches it off or on.
* The dotted marker lines (time and wavelength cursors, the picked scan)
  are only on screen: copied and saved images leave them out.

Tool bar (above the plots)
* Select (default, Esc), Zoom (drag a box; Ctrl + drag does the same with
  any tool; the box stops at the plot edge when the mouse goes past it),
  Pan (drag the view), Full view.
* Background: drag across a chromatogram to set the background range.
* Integration: Auto (every peak in the visible window of the trace chosen
  in the panel), Drag (drag from the start to the end of a peak), Click
  (click near the top of a peak; its limits are found automatically),
  Split (click inside a peak to divide it with a drop line), Delete (click
  a peak to remove it), and the bin icon to clear all peaks. The tools work
  on any chromatogram tile, including mass chromatograms and the PDA
  chromatogram; peaks of several traces can be kept at the same time.
* MS spectra: Mass chrom. (click a peak for its mass chromatogram),
  Measure (click two peaks: distance in m/z; for adjacent charge states the
  charges and the mass are shown; for isotope peaks the charge), Label
  (click a peak to pin its m/z label).
* PDA: Label pins the wavelength and absorbance of a band in the UV
  spectrum.

Peak table
* Trace, retention time, start, end, height, area (signal x seconds, as in
  LabSolutions) and area % within each trace. Selecting a peak shows its
  averaged mass spectra (MS) or its UV spectrum (PDA); Delete removes it.
  Export as .csv.

Every plot: Ctrl + drag to zoom, double click for the full view, right click
for Copy image, Save image as (PNG or TIFF at 300 dpi, PDF, SVG) and data
export. After a click on a plot, Ctrl+C copies it and Ctrl+S saves it. A
short note at the bottom of the window confirms the copy or the saved file;
if it fails, a message gives the reason (details in the log file). Values under the mouse are shown at the top right of each plot. The
wheel scrolls the tiles; drag the gap below a tile to resize it (short tiles
hide their axis titles; saved images keep them).
Images for publications: copied and saved images are drawn at a fixed size,
whatever the size of the tile: 8.5 x 6 cm (one journal column), text 9 pt
(Arial; tick labels about 8 pt, peak labels about 7 pt), 300 dpi. Word and
PowerPoint paste them at exactly that size, so the text stays readable.
Right click > "Image size for copy and save..." changes width, height, text
size and resolution (e.g. 17.5 cm for two columns), or "As on screen". The
margins are fitted so that no axis title or peak label is cut off. PDF and
SVG keep the text as text (editable in Illustrator or Inkscape).
Integration without the panel: right click a chromatogram > "Integrate the
peaks in the visible window..." (trace, minimum height and width), "Export
the peak table...", "Remove every integrated peak"; the peak table has the
same in its right click menu. PDA: right click the map or the chromatogram
> "Map and chromatogram settings..." (bandwidth, colour scale, more
wavelengths, max plot, smoothing, MS overlay and delay) and "Go to a time
and wavelength...".

Checked against LabSolutions exports of KUA525-POZi_027.lcd:
* averaged positive spectrum 2 to 10 min: same m/z values; intensities agree
  (r = 0.9995) when compared peak by peak. LabSolutions exports one point
  per peak, LCMS Analysis keeps the full profile (better for deconvolution).
* PDA chromatogram at 280 nm, 4 nm: r = 0.999 against the LabSolutions plot;
  UV spectrum at 9.398 min: maxima at 221 and 262 nm, within a few mAU.
Not included (LabSolutions functions outside data review): instrument
control, method editing, calibration curves and quantitation, library
search, peak purity, audit trail.

MS reading uses OpenSZRaw 0.2.0 (Apache-2.0, bundled). For LCMS-2020 files
it reports m/z values twice too high; LCMS Analysis detects this from the scan
range stored in the file and corrects it. The PDA decoder and everything
else are part of this package (_portable\unilcms.py, lcms_data.py,
lcms_pda.py, lcms_integrate.py).


LCMS ANALYSIS: FILES OF OTHER VENDORS
-------------------------------------
Open raw data... (Ctrl+O), Open a .D or .raw folder..., or drag the file or
folder onto the window. A file may hold MS data, PDA (UV/Vis) data or both.
* Agilent ChemStation / OpenLab .D folders: DAD spectra (*.uv), single
  wavelength signals (*.ch, shown together as a PDA map when there are no
  spectra), single quad MS (MSD1.MS ...: one scan event per signal, polarity
  from the acquisition method). OpenLab .dx files too.
* Waters .raw folders: MS functions (polarity per function) and the PDA
  function (stored in micro AU, shown in mAU); the PDA analog channels when
  there is no PDA function.
* Thermo .raw files: MS1 scans (the centroids of profile FT scans), events
  by polarity and analyzer; PDA or UV channels.
* mzML and mzXML (MS only).
* ANDI / AIA netCDF (.cdf), exported by Empower, Chromeleon, OpenLab,
  LabSolutions and others: MS data or one chromatogram channel.
* Sciex .wiff / .wiff2 cannot be read: convert them to mzML (ProteoWizard
  MSConvert).
Spectra are averaged on m/z bins (bin width setting) as for .lcd files
(mzML files as in HRMS Analysis).
Other detectors (ELSD, CAD, FID, analog signals) and MRM data are not shown;
the status bar names them.
Agilent and Waters files are read with rainbow 1.5.3 (LGPL 3.0, unmodified,
_portable\wheels; see LICENSES\rainbow_LICENSE_LGPL-3.0.txt), Thermo files
with Thermo's RawFileReader (shipped with the UniDec package, licence in
LICENSES\Thermo_RawFileReader_License.doc), mzXML with pyteomics, netCDF
with SciPy.


COMPARE VIEW (LCMS ANALYSIS)
----------------------------
Compare, the third view in the top bar of LCMS Analysis, puts the
chromatograms of every open file on top of each other. Open the runs
(several at once: Ctrl or Shift + click in Open, or drag them onto the
window), then choose Compare.
* Trace: PDA chromatogram at a wavelength (with bandwidth; or each file at
  the wavelength set in its own PDA view), PDA max plot, TIC, base peak
  chromatogram or a mass chromatogram (m/z with a window in Da or ppm), of
  the positive or negative scans. The first wavelength is the one of the
  PDA view of the file shown; "lambda max of the main peak" sets the
  absorption maximum (above 210 nm) of the tallest peak of the first run,
  with the spectrum just before the peak subtracted.
* Files: every open file is listed in the side panel. Tick the ones to
  compare (Delete or Space unticks the file chosen); Up, Down, Reverse and
  Sort (name, acquisition time, sample name) set the order, the first file
  is drawn at the top. Each file can get its own label, colour and time
  shift; a click on its trace chooses it. Files opened or closed while
  comparing appear in the list or leave it. Runs of the same sample are
  named with their file name too. A blank run (solvent injection) can be
  subtracted from every other trace; it is not drawn itself.
* Processing, the same for every trace: smoothing, baseline (offset,
  straight line between the ends, or rolling minimum), a time window (or
  right click > Use the time range shown).
* Alignment: right click a peak > Align every trace on its peak near here
  (or the Align tool, or type the time): each trace is shifted so that its
  peak near that time lines up. A trace without a real peak there stays
  where it is, and a note says so.
* Layout: stacked (spacing in % of the tallest trace, skew for a waterfall,
  first file at the top or at the bottom), offset (a 2D waterfall as in
  papers: each run a little higher and to the right, overlapping, with a
  real y axis; the front run is drawn over the ones behind it, the time
  axis is that of the runs themselves and the runs behind are cut at its
  end) or overlay. Scale: the same for all (absolute), each trace to its
  tallest peak, or each to its peak at a reference time (right click the
  peak). Y axis: automatic gives stacked traces a scale bar instead of a y
  axis (the axis reads true for the lowest trace only), and the axis with
  its values to offset, overlay and a single trace; or choose scale bar,
  axis or none. On an offset plot, too, the axis values are true for the
  front run only.
* Labels: names at the right or left end of each trace, right of the frame
  at the level of the end of each trace (in its colour; the default of the
  offset layout), a legend, or none (sample name, file name or both). Retention times of the main peak or of
  every peak above a threshold, always the run's own time (also when the
  traces are aligned or shifted). Labels do not sit on top of each other;
  one that finds no room is left out.
* Guide lines: dashed lines at chosen times (Guide line tool, right click,
  or type them); the table gives the area % of the peak at each line, e.g.
  starting material and product in a reaction series.
* Table: per trace the main peak (its retention time in the run), height
  (mAU or counts), area %, the area % at each guide line, time shift, file,
  sample and acquisition time. Notes say when runs differ (method, run
  length) or a file was left out (no PDA data, wavelength outside its
  range, no scans of that polarity).
* Export: the traces as processed (time and value per run) and the table
  as CSV, e.g. for Origin; copy or save the plot (PNG, TIFF, PDF, SVG) at
  journal size, with names, retention times, legend and scale bar placed
  again for the image; Report... makes the comparison report.
* Graph properties (right click the plot): font, text sizes and bold, axis
  titles, x and y range (the x range is then the full view: double click
  goes back to it), tick steps (a step giving more than 50 ticks falls back
  to automatic), ticks inside or outside, minor ticks, frame, grid, colour,
  width and style of each trace, fill opacity, legend position and frame,
  a title above the plot and a panel letter. Apply, OK, Cancel (undoes an
  Apply; also Esc), Reset all. They belong to this plot for the session
  and are not saved; copied and saved images and the Comparison report
  follow them.
* Mass spectra of a peak (m/z in the tool bar, on or off): while it is on,
  a row of two spectrum tiles (ESI+ left, ESI- right) is shown below the
  plot, and a click on a trace (Select tool) picks its peak there: the
  peak top near the click and the start and end of the peak, on the trace
  as processed. The scans of the peak are averaged for the positive and
  the negative ions of that run (or, Spectra in the side panel, only the
  scan nearest to the peak top; the spectrum just before the peak can be
  subtracted), with the bin width of the Mass spectrometry view. Times
  are those of the run itself (without skew, alignment or shift); on a
  PDA trace the MS time is the PDA time plus the MS detector delay of the
  file (PDA view), and the tiles say so. A click on the baseline takes
  the three scans around that time. The strongest ions are labelled with
  their m/z, and a suggested neutral mass is given from the ions of both
  spectra (e.g. [M+H]+, [M+Na]+ and [M-H]-): a suggestion to check, not
  a proof. How it is found (C++ library): the
  peaks of both spectra above 4 % of their base peak, without 13C
  isotopes; each read as [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+, [M+2H]2+,
  [2M+H]+, [2M+Na]+, [M-H]-, [M+Cl]-, [M+HCOO]-, [M+CH3COO]-, [M-2H]2- or
  [2M-H]-; masses that agree within 0.3 form a group of at least two
  ions, scored by their intensities (more when both polarities agree).
  Mobile phase acids (TFA, formic, acetic, difluoroacetic acid) are named
  as such and never taken for the compound; groups that only read the
  same ions again (e.g. as 2M or [M+K]+) are left out; an M+2 partner is
  noted with its intensity ratio and what it suggests (about 1:1 one Br,
  3:1 one Cl, 3:2 two Cl). The line also names the base peaks the
  suggestion does not explain, says "Weak suggestion" when its ions are
  only a few % of the base peaks, and "if ... is [M+H]+" when only one
  ion is found. Subtracting the background (on by default) keeps the
  mobile phase out of the spectra. On a unit resolution instrument the
  mass is good to about 0.3. The peak is marked on its trace (a
  ring at its top, its range shaded), also in copied images and the
  report. Clicking another peak replaces the pick; m/z off removes the
  tiles and the mark. The Comparison report then offers the two spectra
  with a caption (run, time range, delay, averaging, suggested mass).
* Time ranges with m/z on, as on the chromatograms of the Mass
  spectrometry view: drag across a trace (Select tool) to average that
  time range of its run (the trace the drag starts on; times of the run
  itself, PDA traces plus the MS detector delay; a range past the end of
  the run stops there, a range without a scan takes the scan nearest to
  it). Shift + drag sets the background range of that run: it is
  subtracted from the spectra of that run (ranges, single scans and,
  while "Subtract the spectrum just before the peak" is ticked, peaks
  clicked, in place of the spectrum just before the peak). Alt + click
  takes the scan nearest to that time, as a click in the Mass
  spectrometry view. Ctrl + drag zooms a box; Zoom and Pan as always.
  The averaged range is shaded on its trace (a band in the colour of the
  trace, the area under it darker), each background range as a grey
  hatched band on its own trace; the tiles say e.g. "averaged 9.20 to
  9.60 min (range), background 8.80 to 9.00 min", and so does the report.
  Right click the plot: remove the background range of that run, or
  every one; m/z off removes them.
* The spectrum tiles work as the spectra of the Mass spectrometry view
  (the same code): the m/z under the pointer, drag to zoom the m/z axis,
  Ctrl + drag or Zoom for a box, Ctrl + wheel, Pan, double click for the
  full view. While m/z is on the tool bar shows Mass chrom., Measure and
  Label: click a peak in a tile to follow that ion in every run (the
  trace becomes its mass chromatogram, window 0.5), to measure two peaks
  (distance, isotope spacing or charges and mass, for the polarity of the
  tile) or to pin its label (pinned labels stay when another peak of the
  same run is picked, not on the spectra of another run; a new pick
  removes the measurement). The plot itself keeps its
  own tools. Right click a spectrum: mass chromatogram of that ion in
  every run, Label m/z (or remove it), remove every label, remove the
  measurement, Display (profile or sticks; at first as in the Mass
  spectrometry view of the file), show it in the Mass spectrometry view
  of that file (that time range and background), export it (text,
  JCAMP-DX or CSV with the run and times), copy the data, open it in the
  Deconvolute window (negative ions as JCAMP-DX, so it starts in negative
  mode), full view, copy or save the image.
* The side panel of the Compare view is shown or hidden on its own (Show
  panel or F9), apart from the other views.

Region areas and kinetics
* Area tool: drag over a peak (or type the start and end) to integrate the
  same retention time window in every run: straight baseline between the
  ends, positive area only, in signal x s. The table keeps the main peaks and
  guide line area % and adds the region peak, height, area, relative area,
  the 100 % reference run and X. A run without usable data shows "no data".
* X values: one per run in plotted order, separated by commas, semicolons or
  new lines; they stay with their runs when the list is reordered.
* Plot areas vs X, export them to Excel, or fit them (Kinetic fitting:
  straight line, exponential decay, rise to a final level; values with
  standard errors and the half life).


SAVED ANALYSES (LCMS ANALYSIS AND HRMS ANALYSIS)
------------------------------------------------
* After an edit, the analysis is saved by itself in <file>_analysis\
  session.msanalysis next to each edited data file (and on closing). Files
  only opened or zoomed are not written to. Save (top bar) or Ctrl+Shift+S
  saves at once.
* Opening a data file that has a saved analysis asks whether to resume it.
  Open > Open analysis folder / Open analysis project opens one directly, or
  give a .msanalysis file to MSpektra.exe.
* The file keeps the views, results, Compare view, polymer and kinetics
  windows; the raw data stays where it is and is read again. The previous
  save is kept as session.msanalysis.previous.
* Opening an analysis does not change your settings; mass shifts it uses
  that are missing from your list are added to it.
* A second window does not save over an analysis open in another window.


METHOD PRESETS (LCMS ANALYSIS AND HRMS ANALYSIS)
------------------------------------------------
A method is a named set of every setting of a window: the same processing
for every run of a series in one click, and the same settings on another
PC. Method in the top bar (right of the File, MS and PDA or Calibration
chips) shows the method of the file shown ("none" until one is saved or
applied); a click opens its menu:
* Apply to this file, Apply to every open file: the method chosen. The
  controls show its values and everything is computed again as after a
  change in the side panel: chromatograms, averaged spectra (new bin width,
  background on or off), the PDA traces and the UV spectrum, the Compare
  view. The status bar says what was applied, e.g. Method "Gradient A"
  applied to run_B.lcd (Mass spectrometry view, PDA view, Deconvolution,
  Mass shifts, Compare view, MS and PDA link). The settings of the window
  (deconvolution, mass shifts, Compare view, calibration) belong to every
  file: they change with either choice. Integrated peaks and
  deconvolution results are not computed again (Integrate and Deconvolute
  do that, as before).
* Save the current settings as...: a new method from the file shown (a
  method of the same name can be replaced).
* Update "name" with the current settings: the method of the file shown.
  The first line of the menu says "Method: name (settings changed since)"
  when the settings are no longer those of the method.
* Rename, Delete: the open files keep their settings.
* Default when a file is opened: a method applied to every file opened in
  that window (LCMS Analysis and HRMS Analysis each have their own), or None:
  a file opens with the settings of the window, as before. The settings of
  the window (deconvolution, Compare view, calibration, MS and PDA link)
  come with the first file of a window only, so that changes made while
  comparing are kept when more files are opened. With a default that
  calibrates automatically on opening, the first file is calibrated with
  the settings of the method.
* Export to a file, Import from a file...: a method as a small text file
  (name.msmethod.json) to keep with the data or to pass on. Import checks
  every value: a value of the wrong kind or out of range is left out with a
  message (that setting then stays as it is when the method is applied),
  settings this version does not know are ignored, and a method of the
  other window or a file that is not a method is refused with the reason.
  When a method of the same name exists: replace it or keep both.
* What a method contains... lists the settings of the window.

What a method contains
LCMS Analysis:
* Mass spectrometry view: trace (TIC, BPC, mass chromatograms only),
  default mass chromatogram window and its unit (m/z or ppm), mass
  chromatograms stacked or overlaid, smoothing, scale each trace to 100 %,
  link the time axes, subtract the background, bin width, spectrum display,
  peak integration (minimum height and width).
* PDA view: wavelength of the slider (the nearest one recorded), bandwidth,
  colour scale, more wavelengths, max plot, smoothing, scale each trace to
  100 %, overlay of the MS trace, MS detector delay, subtract the
  background, peak integration.
* Link MS and PDA times.
* Deconvolution: every setting of the deconvolution window (method,
  preset, charge and mass ranges, m/z range, minimum intensity, charge
  carrier, baseline, peak shape and isotopes, the settings of UniDec,
  maximum entropy and IsoDec, peaks in the result).
* Compare view: signal, wavelength and bandwidth, each run at its own
  wavelength, polarity, m/z and window of a mass chromatogram, smoothing,
  baseline, time window, alignment time and search window, layout, scale
  and reference time, spacing and skew (stacked and offset), order, colours,
  line width, fill, names of the runs and their text, retention times and
  their threshold, guide lines, y axis, peak table threshold, the spectra of
  the m/z tool (averaging, background, ions labelled).
* Mass shifts (the deconvolution results): the shift list, the tags, at
  most k tags and the other shifts with them, the tolerances and the share
  of the tallest mass used (not Show every matched pair).
HRMS Analysis:
* Mass spectrometry view: as in LCMS Analysis, without the bin width (HRMS
  spectra keep the m/z axis of the instrument).
* Exact mass: ion and tolerance (the formula finder uses both). Formula
  finder: element limits, ranking by the isotope pattern.
* Calibration: reference list (and the custom list), model, HPC order,
  search window, minimum intensity, calibrate automatically on opening.
* Deconvolution: every setting of the deconvolution window.
* Mass shifts: as in LCMS Analysis.
Not in a method: the data of a file (averaged and background ranges, the
m/z lists of mass chromatograms, integrated peaks, labels, the formula
checked, the calibration made on a file, results), the names, colours,
shifts and blank of the runs in the Compare view, the zoom, tile heights,
the side panel, the Graph properties (they last for the session only), the
image size for copy and save, the report settings, and the charge state
tiles of the deconvolution results (shared by both windows). The
Deconvolute window of the start screen is UniDec's own window and keeps
its own settings.
Methods are kept in config\portable_settings.json (key method_presets);
settings files of older versions work as before.


KEYBOARD SHORTCUTS (LCMS ANALYSIS AND HRMS ANALYSIS)
----------------------------------------------------
F1 opens the list of every shortcut (search field; this window or every
window); the help (?) in the top bar still shows the overview.

Anywhere in the window (also while the cursor is in a field)
  Ctrl+O                      open a data file (.lcd; HRMS: a .d folder;
                              .raw and mzML from the Open menu)
  Ctrl+Shift+S                save the analysis
  Ctrl+W                      close the file shown
  Ctrl+R                      create a report
  Ctrl+Tab, Ctrl+Page Down    next open file (Ctrl+Shift+Tab, Ctrl+Page Up:
                              previous)
  Ctrl+1 / Ctrl+2 / Ctrl+3    Mass spectrometry / PDA / Compare view
  Ctrl+E                      export the data of the plot clicked last
                              (chromatogram, spectrum, UV spectrum, the
                              traces of the Compare view)
  Ctrl+D                      deconvolute the spectrum shown (settings)
  Ctrl+I                      integrate the peaks in the visible window
                              (Auto, last settings; MS and PDA views)
  Ctrl+Z                      undo (in a text field: the typing there)
  Ctrl+Y, Ctrl+Shift+Z        redo
  F1                          list of the shortcuts
  F9                          show or hide the side panel

In a view (the focus on a plot or the view, not in a field or a list)
  Esc                         back to the Select tool; cancels a
                              measurement waiting for its second peak
  Home                        full view of every plot
  + / -                       zoom in or out on the time or m/z axis of
                              the plot clicked (the keys that type + and
                              - on your keyboard, e.g. Shift+1 on a Swiss
                              keyboard, or the numeric keypad)
  Left / Right                Mass spectrometry view: the spectra of the
                              previous or next scan (Shift: 10 scans)
  Left / Right, Up / Down     PDA view: move the time or the wavelength
                              slider (Shift: 10 steps)
  Delete                      delete the integrated peak chosen
  Ctrl+C / Ctrl+S             copy the plot clicked as an image / save it

Tools: a single letter while a plot has the focus (click it once)
  S Select, Z Zoom, P Pan, B Background, I Click (integrate a peak near
  its top), D Drag (integrate from start to end), K Split, R Delete (click
  a peak), X Mass chrom., U Measure, L Label, F Formula (HRMS Analysis),
  A Align, G Guide line, M m/z on or off (Compare view; X, U and L while
  m/z is on). The tooltips of the tool bar show the letters.

Typing in a field (a number, a formula, a list of m/z) never switches a
tool, and Ctrl+Z, Ctrl+Y, Ctrl+C, Ctrl+V, Home and the arrow keys keep
their usual meaning there. Lists keep their own keys: the peak table
(Delete), the files of the Compare view (Space, Delete) and the list of
open files (Up, Down, Enter, Space, Delete, Left, Right). Dialogs: Esc
cancels, Enter confirms. The mouse is unchanged (drag, Shift + drag,
Ctrl + drag, double click, wheel, Ctrl + wheel, Alt + click, right click).

UNDO (LCMS ANALYSIS AND HRMS ANALYSIS)
--------------------------------------
Each window keeps its own list of steps (up to 100): Ctrl+Z or Undo in
the top bar goes back one step, Ctrl+Y (or Ctrl+Shift+Z) or Redo goes
forward again. The tooltips of the two buttons and the right click menus
name the step ("Undo: add mass chromatogram m/z 803.20 (+)"); after Undo
the status bar says what was undone, and the file and the view of the
step are shown. A new edit after an Undo removes the steps that could be
redone.

What is undone
* Integration: automatic integration (Auto, Ctrl+I, the settings window),
  peaks added (Click, Drag, Add range), split, deleted, Clear, and the
  settings (trace, minimum height and width), in the MS and PDA views.
* Mass spectrometry view: averaged and background ranges, the time
  clicked (the spectra of one scan), mass chromatograms added, edited or
  removed, labels pinned on the spectra and measurements, and every
  setting of the side panel (trace, mass chromatograms shown, smoothing,
  scaling, linked axes, default window, background subtraction, bin
  width, spectrum display, spectrum to use).
* PDA view: the time and wavelength sliders (a drag or a series of arrow
  keys is one step), averaged and background ranges, more wavelengths,
  max plot, bandwidth, smoothing, colour scale, MS overlay, MS trace, the
  MS detector delay, labels of the UV spectrum.
* Compare view: every setting of the side panel, files ticked or
  unticked, their order, labels, colours and time shifts, the blank, guide
  lines, alignment, time window, reference peak, Graph properties (one
  step per visit of the window), the m/z tool (on or off, the peak or
  range picked, the background ranges, its settings, labels and
  measurements on its spectra).
* HRMS Analysis: the calibration (calibrated, removed or loaded), the
  formula check and its isotope pattern, formula, ion and tolerance, and
  everything of the Mass spectrometry view above.
* Deconvolution: a deconvolution is one step that adds its result; Undo
  removes the result and goes back to the settings used before, Redo
  brings the same result back without computing it again. Results closed,
  hidden or shown, masses added, removed or labelled and "Label all peaks
  in view" are steps too. A setting changed in the settings window becomes
  part of the step of the deconvolution it starts. The deconvolution
  settings are shared by every file and window: undoing another edit of
  the results (a mass, a label) leaves them as they are now.
* "Link MS and PDA times" (right click a chromatogram).
* Method presets: applying a method (to the file shown or to every open
  file) is one step, "apply method name"; Undo goes back to every setting
  it changed, the deconvolution settings and the mass shift settings
  included, and the Method chip shows the method before. A default
  method applied while a file opens is the starting point of that file,
  not a step.
* Mass shifts: on or off for a result, the unmodified species, and the
  settings of the Mass shifts window (one step when it closes). The
  settings are shared by every window: a change made in another window
  is a step of that window only.
Restoring a step sets the controls and computes again what depends on
them (chromatograms, averaged spectra, UV spectrum, traces of the Compare
view), as if the change had been made by hand; a deconvolution is never
computed again (its result is kept with the step).

Steps that are merged: quick repeats of the same kind within 1.5 s (the
arrows of a spin control, the arrow keys moving a slider or the scan, a
slider dragged) are one step. A value typed in a field is one step when
Enter is pressed or the field is left. A dialog (for example Graph
properties or the calibration window) is one step when it closes.

Not undone: zoom, pan, full view and the other changes of the view (tile
heights, closed tiles, the side panel, the view or file shown), opening
and closing files, exports, copied or saved images, and reports. When a
file is closed its steps are removed. Undo and Redo wait while a spectrum
is being averaged (progress window) and, for a step of that file, while
its deconvolution runs: the status bar says so; press again when it is
done (or after Cancel). The Deconvolute window (UniDec's own program) has
no undo.


REPORTS (LCMS ANALYSIS AND HRMS ANALYSIS)
-----------------------------------------
Report... in the top bar (or Ctrl+R, or right click any plot > Create a
report...) asks which report to make; reports that need something not done
yet are greyed with the reason:
  HRMS compound report   chromatograms with the ranges used, spectrum,
                         exact mass (calcd, found, error in ppm, isotope
                         match; from Check a formula), deconvolution
                         masses, calibration.
  LC-MS purity report    PDA chromatogram with the integrated peaks and
                         area %, UV, ESI(+) and ESI(-) spectra of the main
                         peak (the largest), m/z table, TICs. The peaks are
                         those integrated in LCMS Analysis; if none, the
                         report integrates automatically (and says so).
  Deconvolution report   spectrum with fit and charge states, zero charge
                         spectrum and charge states, masses (isotope table
                         for resolved isotopes) compared with the expected
                         masses you type (one per line: "18363 beta-Lg A"),
                         all settings.
  Comparison report      the Compare view: plot as shown, the spectra of
                         the peak picked with m/z, runs, peaks of each
                         trace, region areas, and a methods text.
The HRMS and deconvolution reports are plain: figures, tables and settings,
without explanatory captions.
  SI page                figure a to d (LC-MS: chromatogram, UV, ESI(+),
                         ESI(-); HRMS: chromatogram, spectrum, isotope
                         pattern, zero charge spectrum), caption and a
                         characterisation text to paste into the SI.
The figures are the tiles of the window as you see them (zoom, traces,
labels, markers, fit, charge states, shaded ranges, mass chromatograms),
drawn for the page; a figure without a tile (isotope pattern, calibration
errors, peaks integrated only for the report, a spectrum not shown) is drawn
from the data. The HRMS report holds every result shown in the window
(hidden ones stay out); the charge states only when their tile is shown.
Details: compound, operator, notes (and assignment, expected masses or SI
figure number where they apply); empty fields are left out. The sample
information of the file (sample name, acquisition date, method, vial) is
added by itself. Output: PDF, or Word (.docx) to edit the text. Saved in
the _analysis folder of the data file and opened.
The report libraries (reportlab, python-docx) are in _portable\wheels and
are unpacked into _portable\pylibs the first time (or into your local
AppData if the program folder is read-only). Nothing is installed.


DECONVOLUTION (LCMS ANALYSIS AND HRMS ANALYSIS, IN THE SAME VIEW)
-----------------------------------------------------------------
Right click the spectrum > "Deconvolute this spectrum" (or Deconvolute in
the top bar) opens the deconvolution window. It holds every setting, in
sections that show only what applies to the chosen method:
  Spectrum and data   scan event, polarity (set from the spectrum), charge
                      carrier (positive ions: H+, Na+, K+, NH4+; negative
                      ions: loss of H+, Cl-, formate, acetate; or a custom
                      mass added per charge, with its sign), m/z range,
                      baseline subtraction and its window.
  Method              UniDec, maximum entropy, IsoDec (HRMS).
  Ranges              preset, charge range, mass range, mass step; "Suggest
                      ranges from this spectrum" sets them, and the m/z
                      range, from the data ("Suggest them every time this
                      window opens" does it automatically).
  Peak shape          resolving power or peak width (0 = measured), peak
                      shape, isotopes (envelope or resolved), "Show at
                      resolving power" (display only, see below).
  UniDec settings     iterations, charge, point and mass smoothing,
                      suppression (beta), m/z to mass (integrate,
                      interpolate, smart), data reduction and bin size, data
                      smoothing, intensity threshold, isotope mode.
  Maximum entropy     rounds, minimum rounds, iterations per round.
  IsoDec settings     model, matching tolerance, minimum isotope peaks,
                      similarity threshold, knock-down rounds, thresholds,
                      window.
  Peaks in the result peak window, threshold, number of labels.
The last settings are kept (separately for LCMS and HRMS); "Method
defaults" resets them. Deconvolute closes the
window and runs on the spectrum shown for that scan event.

While it runs: a small progress window with the elapsed time and Cancel.
The calculation runs in a separate Python process, so the program stays
responsive and Cancel stops it at once (the process is ended; the next run
starts a new one). Settings that would take hours or run out of memory are
refused before anything starts, with the reason: a mass step below 0.001 Da,
or a mass x charge grid that would need more than half of this computer's
memory (the message says how much it would need). Example: 0.0001 Da over
2600 to 2750 Da is 1.5 million points and no instrument resolves it; 0.01 Da
is already 20 points across a peak at resolving power 80 000.

Speed and multicore: the UniDec engine runs in msengine.dll and maximum
entropy in msmaxent.dll (no Python fallback); both use all cores (OpenMP; no installation or admin rights
needed). What made runs slow before was mostly the grid: maximum
entropy now works only on the m/z windows the chosen masses can reach at
the chosen charges, and stops as soon as the fit matches the noise. G4COL
(2.66 kDa, charges 2 to 3, 0.01 Da): about 8 to 12 s for the calculation;
the first run of a session adds the start of the worker process (a few
seconds on Windows).

Results: the fit (orange) and the charge states of the selected mass
(tags just above the plot) on the m/z spectrum; a row with the zero charge
mass spectrum and the charge states; the mass table. Every run adds a new
result row above the older ones (newest at the top, at most 8 kept); the
table and the fit belong to the active result (the newest, or the row
clicked last, marked "shown in the table"). Right click a result > Hide
this result or Hide all results; right click a result plot >
Deconvolute again with other settings, Zoom to the masses found, Full view
(whole mass range; Ctrl + drag zooms at any time, also after a UniDec run).
On the zero charge mass spectrum, right click a peak > Label ... (keep the
label) keeps its mass label (also in the report); right click where no
mass was found > Add the mass here adds the highest point near the click
to the masses (table row marked "added", _peaks.csv, report); right click
it again > Remove the added mass.
* Intensities are absolute: each mass has the summed height of its charge
  states in the m/z spectrum (the convention of Bruker DataAnalysis), and
  the charge state plot shows the height at each charge.
* The area under the zero charge spectrum is lightly shaded; with resolved
  isotopes every isotope peak is labelled with its mass (4 decimals), and
  the plot title gives spectrum, method, m/z range and mass step, as in
  DataAnalysis.
* Why earlier results looked "too smooth": maximum entropy starts flat and
  becomes sharper with each round; the old version stopped too early (the
  empty baseline counted in the fit quality) and used a resolving power
  that was measured too low (61 000 instead of about 84 000 on G4COL, which
  also corrects the figure in the previous README). Both are fixed. "Show
  at resolving power" broadens the display on purpose if you prefer a
  smoother look (e.g. 10 000); heights are kept.
* Methods:
    UniDec (Bayesian)    the UniDec engine with all its settings; with a
                         mass step below 0.1 Da the isotope peaks are
                         listed as well.
    Maximum entropy      charge state deconvolution by maximum entropy with
                         one charge envelope shared by all masses, on a log
                         m/z grid (TOF: peak width grows with m/z), written
                         for MSpektra from the published principle.
    IsoDec               monoisotopic masses and charges from resolved
                         isotope patterns (HRMS spectra only).
  A better algorithm? For charge deconvolution of ESI spectra, UniDec and
  maximum entropy are the standard methods; the maximum entropy here now
  reproduces DataAnalysis (below). FLASHDeconv (OpenMS) is strong for many
  overlapping proteoforms but needs a large extra installation (pyopenms);
  it is not included.
* Comparison with Bruker DataAnalysis MaxEnt, same spectrum (G4COL,
  maXis II, 6.90 to 7.00 min, charges 2 to 3, 0.01 Da, calibrated):
      DataAnalysis            MSpektra, maximum entropy, resolved
      2663.2901  6.80e6       2663.2935  5.68e6
      2664.2960  8.90e6       2664.2962  7.86e6
      2665.2964  6.85e6       2665.2992  6.64e6
      2666.3001  3.95e6       2666.3014  3.13e6
      2667.3018  1.80e6       2667.3041  1.25e6
  Masses agree within 1.3 ppm, the isotope patterns with R2 0.98; the
  absolute heights are 3 to 30 % lower here, most for the weak isotopes.
  The 2+ ion carries an interfering peak at m/z 1333.66 in the data, which
  affects the M+2 isotope in both.
* Other tests (the true masses are known): beta-lactoglobulin (single
  quadrupole): UniDec 18280 and 18366 Da, maximum entropy 18282 and 18367
  Da (variants B and A). Synthetic isotope resolved spectrum (monoisotopic
  18250.000 and 18330.000 Da): IsoDec 18250.002 and 18329.984 Da.
* Presets set the charge and mass ranges (peptides, denatured proteins,
  antibodies, native complexes, small molecules); every field can be
  changed.
* Every run is saved in the output folder of the data file:
  <spectrum>_<method>_mass.txt, <spectrum>_<method>_peaks.csv and the input
  spectrum (UniDec also writes its usual _unidecfiles folder).


MASS SHIFTS (DECONVOLUTION RESULTS)
-----------------------------------
Which deconvoluted masses are the protein, the protein with a tag, with a
tag and water, with Na, oxidised, carbamylated, and how much of it is
conjugated. Right click the zero charge mass spectrum of a result > Mass
shifts (a tick: on or off for that result), or Mass shifts... in the
Deconvolution section of the side panel (the settings window, for the
result shown in the table).

On the plot: a bracket joins two masses whose difference matches a shift
of the list (green) or one or more tags with up to two other shifts
(violet), labelled with the name and the difference ("Na +22.0", "PDI +
water +861.0"); a mass that nothing explains gets a red label above its
peak with its difference from the reference ("+429.1 unknown"; the 8
tallest in view, all of them in the table). Dotted lines lead down to the
peaks. By default each
mass is joined once: to the reference, or to its neighbouring species
(the nearest mass it follows from, e.g. tag + water from the tag); right
click > Show every matched pair draws every explained pair (40 at most).
Brackets sit just above what lies below them (peaks, labels, other
brackets); when they need more room the plot makes room above the peaks;
the ones that still do not fit are counted ("2 more in the table"), tags
before the other shifts. They follow the zoom and are part of copied and
saved images.

Next to the plot (one tile, right click > Close this tile turns the
mass shifts off for the result):
  Mass shifts              one row per pair: from, to, difference, match
                           (the other match within the tolerance in
                           brackets, e.g. "sulfate adduct (or phosphate
                           adduct)"; one of the same mass, such as the
                           loss of water read backwards for water, is the
                           same match and not shown), error, and what the
                           To mass is ("reference + PDI + water"). The
                           pairs with the reference, the links drawn, and
                           with Show every matched pair all of them.
                           Click a row: its To mass is selected in the
                           plot and the mass table. Right click a row >
                           Use ... Da as the unmodified species.
  Degree of conjugation    (with a tag) the species with 0, 1, 2, ... tags
                           (each with its matched shifts: tag + water and
                           tag + Na count as one tag), their share by
                           height and by area, the masses counted, and the
                           tags per molecule (by height, by area). Unknown
                           masses are left out (the line above says how
                           many).

The reference ("unmodified species") is chosen by itself: the tallest
mass, except with a tag when the tallest mass is a conjugate (common at a
high degree of labelling): then the lightest species of the tag series
below it (one, two, ... tags fewer than the tallest; of several, the one
with the most tags fewer, then the fewest other shifts, then the tallest)
becomes the reference, so that the degree of conjugation counts the
unmodified protein. Example: 18363, 19206 and 20049 Da (heights 40, 100,
30) with a tag of 843 Da: reference 18363 Da, 0.94 tags per molecule
(with the tallest as the reference the unmodified protein was left out:
0.23). The table says how the reference was chosen ("reference 18363.0 Da
(the lightest of the tag series)", "(tallest)" or "(chosen)"). A fully
labelled sample (only one and two tags present) keeps the tallest as the
reference: nothing shows that it carries a tag. Right click a peak > Use
18385.0 Da as the unmodified species sets it by hand (also in the settings
window and the table); Automatic reference goes back. Every mass gets its
composition relative to the reference, directly or through other masses
(e.g. with the conjugate as the reference: the protein is "reference -
PDI", its oxidised form "reference - PDI + oxidation").

The Mass shifts window (kept in the settings):
  Tolerance       envelopes (average masses): empty = 1.5 Da, at least 1.5
                  mass steps; isotope resolved results: empty = 10 ppm. In
                  Da or ppm (of the heavier mass of a pair).
  Masses used     masses below this share of the tallest one are left out
                  (default 1 %; 0 = every mass of the table).
  Unmodified species  automatic (see above) or one of the masses of the
                  result.
  Tags            up to three, each with a name, its average mass and its
                  monoisotopic mass (for isotope resolved results; empty:
                  the average one); at most k tags per molecule (default 4)
                  with up to 0, 1 or 2 other shifts.
  Shifts          the list (name, difference of the average mass and of
                  the monoisotopic mass, ticked: in use); Add..., Edit...
                  (or double click), Remove, Defaults. Differences are
                  signed (a loss is negative); a difference also matches a
                  shift the other way round ("-Na": the reference carries
                  the Na). Add your own, e.g. a variant of the protein
                  (beta-lactoglobulin B is A - 86.09 Da) or a buffer
                  adduct.
  Display         Show every matched pair.

Default shifts (average / monoisotopic, Da): Na (Na - H) 21.9818 /
21.98194, K (K - H) 38.0904 / 37.95588, oxidation 15.9994 / 15.99491,
double oxidation 31.9988 / 31.98983, water 18.0153 / 18.01056, loss of
water -18.0153 / -18.01056, disulfide (-2 H) -2.0159 / -2.01565,
carbamylation (urea) 43.0247 / 43.00581, acetylation 42.0367 / 42.01056,
methylation 14.0266 / 14.01565, formylation 28.0101 / 27.99491,
phosphorylation 79.9799 / 79.96633, phosphate adduct (H3PO4) 97.9952 /
97.97690, sulfate adduct (H2SO4) 98.0785 / 97.96738, TFA adduct 114.0233 /
113.99286, hexose (glycation) 162.1406 / 162.05282, beta-mercaptoethanol
76.1176 / 75.99829. Each is summed over its elemental composition:
monoisotopic masses from AME 2020 (Wang et al., Chin. Phys. C 2021, 45,
030003), average masses from the IUPAC 2005 atomic weights (Wieser, Pure
Appl. Chem. 2006, 78, 2051), the values Unimod and ExPASy use.

Which masses: proteins (isotope envelopes) are compared by their average
masses. Isotope resolved results by monoisotopic differences: IsoDec gives
monoisotopic masses; resolved species (maximum entropy or UniDec with
resolved isotopes) are compared by their most abundant isotope, and a
match may then be one isotope (1.00335 Da) off, written e.g. "K (isotope
+1)"; when the peaks of a resolved result are its isotope peaks, those 1
or 2 isotopes above another peak are grouped with it first.

How a match is chosen (src/shifts.cpp and ms_shifts.py describe every
step): a difference is explained by one shift, or by k tags of one kind
with up to two other shifts (a shift may repeat; two that cancel, such as
water and loss of water, do not count), within the tolerance. Of several,
the one with the fewest terms wins, where each half tolerance of error
counts as one more term: "PDI + 2 Na" (exact) beats "PDI +
carbamylation" (1 Da off at a tolerance of 1.5 Da), "water" beats a three
term match 0.2 Da nearer. Compositions grow from the reference outwards
through the matched pairs, the simplest first; a mass reached by nothing
is unknown. No composition has more tags of a kind than the maximum or
more than one shift more than a direct match could have. The result is
the same every time.

Reports: the deconvolution report has a section Mass shifts (the table and
the degree of conjugation) when they are shown for the result; the HRMS
compound report adds them after the mass table of each result shown.
Copied and saved images and the report figures show the brackets.

Tested (tests/check_shifts.py in the source package): the library and its
Python reference give the same pairs, compositions and numbers bit for
bit on 600 random series (proteins, peptides, small molecules, isotope
peaks, noise, unknown masses) and on known cases: a beta-lactoglobulin
conjugate series (18363 Da, +843 tag, +861 tag + water, +22 Na, +43
carbamylation, +16 oxidation, +883 tag + water + Na, and an unknown +429.1:
all found, 55.6 % unconjugated and 44.4 % with one tag by height, 0.44
tags per molecule), the conjugate as the reference, losses, two Na, isotope
resolved small molecules at 10 ppm. On real data: UniDec of B Lac_013.jdx
gives 18366.6 Da (variant A, the reference) and 18280.7 Da, "variant B
-85.9" with the variant added to the list (+0.19 Da off); the 12 kDa
protein of Protein POS.d shows its oxidation (+15.9 Da).


HRMS ANALYSIS: DATA FILES
-------------------------
Reads Bruker .d folders with analysis.baf (maXis, maXis II, impact,
compact, micrOTOF, apex/solariX) with Bruker's Baf2Sql library (bundled in
_portable\baf2sql), .d folders in the newer TSF format (analysis.tsf and
analysis.tsf_bin, otofControl 6 and later) with Bruker's TDF SDK (bundled
in _portable\timsdata), timsTOF .d folders (analysis.tdf and
analysis.tdf_bin; MS1 frames with the ion mobility dimension summed, MS/MS
frames left out; same SDK), and mzML files (profile or centroid; e.g.
converted with ProteoWizard msconvert). Open:
the Open button offers "Bruker or Agilent .d folder", "Waters .raw folder"
and "mzML or Thermo .raw file"; or drag the folder or file onto the window.
Other vendors (MS1 scans; MS/MS scans are counted, not shown):
* Agilent MassHunter .D folders (TOF, Q-TOF): AcqData\MSScan.bin with
  MSProfile.bin (profile, run length or LZF compressed) and/or MSPeak.bin
  (centroids), read with the rainbow library (rainbow-api 1.5.3, LGPL,
  _portable\wheels). The m/z values are those of the calibration stored
  with each scan (MSMassCal.bin, or DefaultMassCal.xml). One event per
  polarity (and per collision energy when there are several).
* Waters .raw folders (Synapt, Xevo, also quadrupoles) with Waters'
  MassLynx library (MassLynxRaw.dll in the private Python): the m/z values
  with the calibration of the file; lock mass correction is not applied
  (use the calibration of HRMS Analysis). One event per MS function.
* Thermo .raw files (Orbitrap, LTQ FT) with Thermo's RawFileReader: profile
  scans with their centroid lists; one event per polarity and analyzer
  (FTMS, ITMS) and SIM apart from full scans.
* Sciex .wiff: not readable without Sciex's licensed library; convert to
  mzML with ProteoWizard msconvert.
The spectra of these files are read in Python; averaging, mass
chromatograms and everything after run in msengine as for Bruker data.
The first time a .d folder is opened, Baf2Sql writes a small index file
(analysis.sqlite) into it; for a write-protected folder a temporary copy is
used. The m/z values are those of the last calibration saved with the file
(e.g. in DataAnalysis); the calibration made here comes on top.

One view, as in LCMS Analysis, with high resolution defaults
* TIC, base peak and mass chromatograms per polarity; the default mass
  chromatogram window is +-0.02 m/z (each one can have its own, in m/z or
  ppm). Spectra keep the instrument's profile m/z axis (no binning); labels
  and read-outs show 4 decimals. MS/MS scans are counted but not shown.
* Exact mass: right click the spectrum > "Check a formula": type the
  formula of the neutral molecule and choose the ion ([M+H]+, [M+Na]+, [M+K]+, [M+NH4]+, [M+2H]2+, [M]+, [M-H]-,
  [M+Cl]-, [M+HCOO]-, ...). Check shows the calculated m/z, the measured
  peak with the error in ppm, and the calculated isotope pattern drawn over
  the spectrum, with an isotope match score.
* Formula tool (tool bar) or right click a peak > "Find formulas": elemental
  formulas within the tolerance (element limits editable; C, H, N, O, S, P,
  F, Cl, Br, Na, Si, B), ranked by mass error and isotope pattern; RDB and
  plausibility rules (Kind and Fiehn) filter impossible formulas. Double
  click a result to draw its pattern.

Calibration (the calibrant is in the same run)
* What it does: the calibrant (sodium formate by default, injected at the
  start of the run) gives peaks of exactly known m/z. Their measured m/z are
  compared with the true values, a smooth correction curve (the model) is
  fitted through the differences, and the correction is applied to every
  spectrum and mass chromatogram of the file.
* How:
  1. Drag across the calibrant peak on the chromatogram. When a file opens,
     the calibrant scans are looked for and marked green on the TIC; right
     click the chromatogram > "Show the calibrant spectrum" selects them.
  2. Right click the spectrum > "Calibrate on this spectrum" (or right click
     the chromatogram > "Calibrate on the calibrant").
  3. The calibration window shows at once the calibrant spectrum (green:
     reference ions found and used), the mass error of every reference ion
     as recorded (grey) and calibrated (blue) with the model (orange), the
     table of reference ions (double click one to leave it out), and all
     models compared on the same peaks. Choose calibrant, model and search
     window; press Calibrate: it is applied to the whole file. The window
     has the same scrolling tiles as the main window (journal proportions,
     drag a gap to resize one tile, right click > Close this tile).
     Removing calibrants from the fit: select rows in the table (Ctrl or
     Shift + click for several) and press Delete or right click > "Remove
     the selected calibrants"; or right click a point in the error plot or
     a green peak in the spectrum > "Remove m/z ... from the fit"; or right
     click > "Remove calibrants with an error above..." (ppm). "Use every
     calibrant again" brings them all back. The fit and the models table
     are recalculated at once.
  4. Close: the main window shows the calibrated spectra; the top bar shows
     the result (e.g. "0.54 ppm, HPC"); click it to open the window again.
  If the window is closed with a fit that was not applied, it asks.
* Models: single point (ppm shift), linear, quadratic, cubic, TOF (square
  root) and HPC (default): quadratic, then a polynomial correction of the
  errors it leaves inside the calibrated range, as in DataAnalysis
  "Quadratic + HPC" (Gobom et al., Anal. Chem. 2002, 74, 3915). Bruker does
  not publish its exact HPC algorithm, so the numbers can differ slightly
  from DataAnalysis. The order of the correction is chosen automatically by
  leave-one-out cross-validation (or set 1 to 8); outside the calibrated
  range the correction is held at its edge value.
* Cross-validated error (CV): each reference ion predicted by a fit made
  without it. The ordinary RMS always gets smaller with more terms; the CV
  error shows how well a model does for peaks between the calibrants, so
  use it to choose the model.
* Calibrant lists: sodium formate (positive, negative), ESI-L tune mix
  (positive, negative), caesium iodide, or a custom list (typed or from a
  file). Every calibration is saved as <name>_calibration.txt (report) and
  .json in the output folder; "Load calibration..." applies a saved one to
  another file. "Calibrate automatically when a file is opened" does it on
  the calibrant scans found, with the last settings.
* Tested on the sodium formate spectrum of a maXis II run (G4COL, 20
  reference ions from m/z 363 to 1655): as recorded -2 to -14 ppm (11.85
  ppm RMS). RMS / cross-validated: quadratic 0.84 / 1.17 ppm, cubic 0.50 /
  0.69 ppm, TOF 0.87 / 1.22 ppm, HPC (order 3) 0.54 / 0.94 ppm. On this
  spectrum cubic predicts best (lowest CV); the Models table shows this for
  every file, and the model chosen is kept for the next time.

Details tile: the isotope peaks of the selected mass (m/z, height, area,
charge) and the measured against the calculated isotope pattern.

Polymer analysis (top bar, after a deconvolution): mass distribution (Mn,
Mw, dispersity, by area or height), repeat unit finder, chain lengths and
end groups from the species of the active result.

Not included (Bruker DataAnalysis functions): methods and automation,
compound lists, SmartFormula scoring, MS/MS spectra, mobility data.


DECONVOLUTE WINDOW: ADD-ONS (UniDec's own code is unchanged)
-----------------------------------------------------------
The Deconvolute tile opens the UniDec window (renamed in the interface).

Right-click menu on every plot
* Save plot as... (PNG or TIFF at 300 dpi, PDF, SVG, EPS; white
  background, cropped to the plot). The file name is suggested from the
  data file, e.g. "B Lac_013_mass_distribution.png", in the data folder.
* Copy plot to clipboard, Export plot data as text, Reset zoom, Clear labels.
* The hidden right-click actions of UniDec are now menu items instead of
  firing on every right-click:
    Mass Distribution plot: "Integrate peaks (shade peak areas)"
    MS Data plot:           "Set m/z range to current view (reprocess)"
  Ctrl + right-click still does what it did before (e.g. Ctrl + double
  right-click on the MS Data plot removes the zoomed region).
  The UniChrom chromatogram keeps its right-button zoom.
* Middle-click on a plot still opens UniDec's own save dialog.

Modern interface
* Sharp text and plots on scaled displays (125 %, 150 %): UniDec now tells
  Windows that it handles display scaling itself, so Windows no longer
  stretches (blurs) the whole window. The log in the logs folder shows the
  display scaling that was detected.
* Deconvolute window and MetaUniDec: a toolbar with labelled buttons (Open,
  Process, Deconvolute, Pick peaks, Run all, Replot) and the file name,
  number of points, R2 and UniScore on the right. Section headers are clean
  bars with a thin colour mark; a click anywhere on a header opens or closes
  it. The main buttons (Process Data, Run Deconvolution, Peak Detection) are
  blue. The bottom quick buttons are named after what they open
  (Processing, Deconvolution, Peaks instead of Blue, Yellow, Red).
* Deconvolute window: every plot sits on a white card with a title. The titles are part
  of the window only, not of saved or copied figures.
* Peak list: white rows with a small marker in the peak colour instead of
  fully coloured rows.
* Windows reopen at their last size and position; the Open dialog starts in
  the folder of the last opened file.
* Start screen: three tiles (LCMS Analysis, HRMS Analysis, Deconvolute);
  arrow keys and Enter also work. Right click the background for the
  classic tools.
* Headings, toolbar and launcher use the Plus Jakarta Sans font (bundled in
  _portable\fonts, SIL Open Font License, see LICENSES); it is loaded for
  this program only and not installed on the computer. Input fields keep the
  Windows system font.
* All buttons keep their functions and Ctrl shortcuts. The add-on is
  _portable\unidec_theme.py; without that file the program starts with
  UniDec's original interface (window sizes are stored in config\portable_settings.json).

Modern plot style
* Plots are drawn in a closed box, with inward major and minor ticks on
  the bottom and left axes (where the numbers are). Intensity axes read
  "Relative intensity (%)" with 0-100 ticks (only the tick labels are
  scaled; data, exports and peak values are unchanged). There is some
  head room above the tallest peak.
* Dark grey data line, orange fit, colour-blind friendly peak colours
  (blue, orange, green, ...) with white-edged markers, lighter integration
  shading, frameless legend, italic m/z. Axis numbers no longer switch to
  "+1.83e4" style offsets when zoomed in.
* The peak colours are used only while the peak colour map in UniDec is the
  default ("rainbow"); any other colour map you choose there is respected.
* Right-click any plot > "Classic UniDec plot style" switches back to the
  original look (and again to return). The choice is remembered.
  The 2D plots (m/z grid, mass vs charge) keep their original look.

Movable labels
* Drag a label with the left mouse button to move it. Once it is moved away
  from its peak, it gets a white box and a thin leader line to the peak.
  Moved positions survive redrawing (integration, Ignore, zoom, style
  switch), are kept in saved figures, and are remembered when the same
  peaks are labelled again.
* Double-click a label to edit its text (several lines are allowed, e.g.
  "18,367 Da" and a name below). Right-click a label for "Edit label
  text...", "Reset label position" and "Delete this label".
* Everywhere else the left mouse button zooms as before.

Labels stay until you remove them
* Labels made from the peak-list menu (masses, areas/intensities, names,
  differences, charge states) survive redrawing: integrating, colour
  changes, Repopulate, zooming.
* They are removed by: "Clear labels" (right-click on the plot, or the new
  "Clear Labels" entry at the bottom of the peak-list menu); choosing
  another label type (it replaces the previous set, as before); "Ignore" or
  "Isolate" (only the labels of the hidden peaks go); and new results
  (opening a file, Process Data, Run UniDec, Peak Detection, All).

Protein LC-MS preset
* File > Presets > Custom > "Protein LC-MS 10-30 kDa (denatured).dat":
  m/z 600-2000, mass 10-30 kDa, sample mass every 1 Da, peak detection
  range 20 Da, threshold 0.05, integration window -40 to +40 Da.
  On B Lac_013.jdx the UniDec defaults (10 Da sampling, 500 Da peak range)
  report one protein peak at 18,370 Da plus low-mass artefacts; this preset
  reports 18,367 Da (variant A) and 18,280 Da (variant B, 34 % height,
  35 % area). The 500 Da default peak range merges species closer than
  500 Da (variants, small adducts).
* Note: UniDec 8.2.1's own default integration window is -1 to -1 Da (the
  box to change it was removed in 8.1.3), so "Integrate peaks" gives areas
  of 0 with the default settings. The preset fixes this for its runs.

Faster start
* Python's compiled files are included (no one-time compilation on the
  first start, also on read-only or pendrive copies).
* The Thermo .raw reader (which starts .NET) now loads only when a .raw
  file is opened, instead of at every start.
* The log (logs folder) shows "Libraries loaded in X s" for each start.


JCAMP-DX (.jdx) ADD-ON
----------------------
UniDec in this package also opens JCAMP-DX mass spectra (.jdx, .dx,
.jcamp), e.g. the averaged LC-MS spectra exported by Shimadzu LabSolutions.
* Open them like any other file (File > Open, the folder button, or drag
  the .jdx onto the UniDec window). No separate conversion is needed.
* UniDec saves the imported spectrum as text in
  <name>_unidecfiles\<name>_rawdata.txt next to the .jdx.
* Works in UniDec (also batch mode: drop several files on the window) and
  MetaUniDec (Add Data Files). Polarity (ESI+/-) and the retention-time
  window are read from the file header.
* Convert_JDX_to_TXT.bat: drag one or more .jdx files (or a folder) onto it
  to save <name>_unidec.txt (tab separated m/z and intensity) next to each
  file, same output as the separate JDX converter.
Tested on a Shimadzu JCAMP-DX 4.24 file: the imported data are identical to
the separate converter's text output and give identical deconvolution
results. Compressed JCAMP-DX (X++(Y..Y)) data, peak tables, several spectra
per file and NTUPLES files are also supported.
The add-on lives in _portable\unidec_jcamp.py; UniDec's own code is unchanged.

START-UP SPEED
--------------
The start screen needs only the window library (about 1 to 2 s); LCMS and
HRMS Analysis load NumPy and Matplotlib, the Deconvolute window loads
UniDec's interface (the largest part). These libraries load in the
background while the start screen is shown. Each log file (logs\) lists
the times: "Ready for the first window after", "start screen shown in",
"... libraries loaded in the background in", "window built in", "read
<file> in" and "UniDec interface loaded in".

Loading is limited by reading about 11 000 files, not by the processor
(the deconvolution itself uses every processor core). What makes it slow:
  * OneDrive: in a OneDrive folder, files are fetched from the cloud the
    first time on each computer ("Files On-Demand"), and the sync client
    checks every file that is opened. Keep MSpektra in a local folder
    that is not synchronised, e.g. C:\MSAnalysis, or right click the
    folder > OneDrive > "Always keep on this device". To compare: copy the
    folder to C:\MSAnalysis, start both copies twice and compare the times
    in their log files.
  * The virus scanner checks every program file the first time it is
    opened; the second start is faster. An exclusion for the folder (ask
    your IT) removes most of this.
  * The very first start on a computer also builds Matplotlib's font list
    (config\matplotlib, kept for later starts) and compiles the add-ons.
Do not pack the program into a single .exe: such files unpack everything
to a temporary folder at every start, which is slower.
Settings (config\portable_settings.json): "start_maximized": false opens
the windows at their normal size; "prestart_deconv_worker": false starts
the deconvolution process only at the first deconvolution.

IF SOMETHING GOES WRONG
-----------------------
* Start  MSpektra (console).bat  instead. It shows all messages in a black
  console window; the window stays open after a crash so the error can be
  read (or photographed / copied).
* UniDec opens but nothing can be processed: the program folder or the
  folder holding your data file is write-protected. UniDec writes a
  "<datafile>_unidecfiles" folder next to every data file, so keep data in
  a folder you can write to. (Program-folder protection is handled
  automatically since the launcher update of 24.09.2026.)
* When started with MSpektra.exe or MSpektra.bat, the same messages are written to the
  "logs" folder (last 20 sessions are kept), or to
  %LOCALAPPDATA%\UniDecPortable\logs if the program folder is
  write-protected. Native crashes are recorded there as well.
* "unidec Run Error: 3221225477" (a memory crash in the deconvolution
  engine): the developer's advice is Advanced > Reset To Factory Default, and
  deleting the "<datafile>_unidecfiles" folder next to the data file, then
  rerunning. If it persists, it is data/parameter specific and belongs in
  a report to the developer (GitHub issues), with the data file.
* To reset only this package's plotting/caching settings, delete the
  "config" folder (it is recreated on the next start).
* If an add-on misbehaves, the log says which one; each add-on is a file in
  _portable\ (unidec_ui_addons.py, unidec_plotstyle.py, unidec_theme.py,
  unidec_jcamp.py, unidec_fast.py, unilcms.py, hrms.py, deconv_tab.py) and
  the program starts without it if that file is removed.
* The plot style choice, window positions and the last data folder are
  stored in config\portable_settings.json (delete it to reset them).
* If Windows shows a security prompt for MSpektra.exe/.bat or python, that is the
  "downloaded from the internet" flag. Step 1 above prevents it; the
  package also removes this flag from its own files on the first start.


WHAT IS DIFFERENT FROM THE OFFICIAL ZIP
---------------------------------------
* No PyInstaller-packed GUI_UniDec.exe. The program starts through the
  signed python.exe/pythonw.exe from the Python Software Foundation, which
  security software and locked-down institutional PCs block less often
  (the developer's release notes mention antivirus blocking of the exe).
* The Windows "Mark of the Web" download block is removed automatically on
  first start (this block is what the official instructions work around
  with "Unblock", and what breaks Thermo .raw loading when it is missed).
* Isolated from other Python installations: no PATH, PYTHONPATH, per-user
  site-packages or matplotlib settings from Anaconda etc. can interfere.
* Readable logs of every session.
* The deconvolution engine (unidec.exe / unideclib.dll with Intel MKL, used
  by the UniDec windows of the start screen; LCMS Analysis, HRMS Analysis and
  the Deconvolute window run the engine in msengine.dll),
  the Thermo RawFileReader and the Waters MassLynx DLLs are the developer's
  own binaries, unchanged.


FILE FORMATS
------------
Supported: Thermo .raw, Waters .raw folders, mzML, mzXML, text (.txt/.dat/
.csv: two columns m/z and intensity, or three columns for IM-MS), .npz,
I2MS (.dmt/.i2ms), UniDec HDF5 files.
Deconvolute window, not in UniDec 8.2.1 (removed upstream in 8.1.0):
Agilent .d and Sciex .wiff direct import, Isotope mode. (Bruker .d: use
HRMS Analysis, which can send spectra to the Deconvolute window.) Convert those files to mzML with ProteoWizard
MSConvert first.
Thermo .raw reading uses the .NET Framework 4.x that ships with Windows
10/11 (nothing to install).
PDF report generation additionally needs MiKTeX (pdflatex) on the PATH,
exactly as for the official version.
The optional MassQL query tool is not included.


CONTENTS
--------
MSpektra.exe             start MSpektra (no console window); the only
                            program file in this folder, the shortcuts and
                            the taskbar point to it
MSpektra.bat             same, fallback if MSpektra.exe is blocked
MSpektra (console).bat   start with a console for troubleshooting
Create_desktop_shortcut.bat make desktop and Start menu shortcuts
Convert_JDX_to_TXT.bat      convert .jdx files to text files (drag and drop)
_portable\                  launcher and add-ons (baf2sql\ and timsdata\: Bruker readers)
_portable\msengine\         C++ libraries: msengine.dll (readers, deconvolution,
                            calculations), msmaxent.dll (maximum entropy),
                            mskinetics.dll (kinetic fits), mspolymer.dll
                            (polymer analysis); sources in the repository (cpp\)
python\                     private Python 3.12.10 with all packages
  python\Lib\site-packages\unidec\bin\   UniDec engine, presets, settings
logs\, config\              created on first start
backups\                    earlier copies of replaced files
                            (old_launcher\UniDec.exe: the old launcher, not used)

Main package versions: unidec 8.2.1, numpy 2.4.6, scipy 1.18.0,
pandas 3.0.3, matplotlib 3.11.0, wxPython 4.2.5, numba 0.65.1,
h5py 3.16.0, pythonnet 3.1.0, pymzml 2.6.1, openszraw 0.2.0, olefile 0.47 (dependency versions pinned to
what was current when UniDec 8.2.1 was released, June 2026).


LICENSES
--------
See the LICENSES folder. UniDec is distributed under its own BSD-style
license (redistribution permitted with the copyright notice and the
citation request above). Python is under the PSF license. OpenSZRaw is
under the Apache License 2.0 and olefile under the BSD 2-clause license
(texts in LICENSES). The Thermo
RawFileReader and Waters MassLynx SDK components are under their vendors'
licenses. Bruker's Baf2Sql library (_portable\baf2sql) is Bruker software,
freely distributed by Bruker (copy from the pyBaf2Sql project, Apache-2.0);
it uses the Intel MKL and other components listed in
LICENSES\Bruker_Baf2Sql_THIRD-PARTY-LICENSES.txt. Bruker's TDF SDK library
(_portable\timsdata\timsdata.dll, for analysis.tsf and analysis.tdf) is
distributed under the Bruker Software License Agreement in
_portable\timsdata\LICENCE-BRUKER.txt,
with THIRD-PARTY-LICENSE-README.txt and redist.txt next to it (copy from the
opentims_bruker_bridge package). This software uses TDF Software
Development Kit software. Copyright (c) 2019 by Bruker Daltonik GmbH. All
rights reserved. The other Python packages keep their own license files inside
python\Lib\site-packages\<package>.dist-info.
