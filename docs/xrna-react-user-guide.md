
# XRNA-REACT User Guide

## About XRNA-REACT

XRNA-REACT is a web-based tool for editing and visualizing RNA secondary structures. It is an upgraded version of the original XRNA software developed by Harry Noller (available as a depreciated [Java applet](http://rna.ucsc.edu/rnacenter/xrna/xrna.html)). Please note that XRNA can be used with a trackpad but works best with a mouse, as several functions require a mouse wheel.

## Basic Navigation

1. To begin working with XRNA, you need an input file:

:::{note}
```{eval-rst}
Download an :download:`example XRNA </files/5S_E_coli.xrna>` or an :download:`example JSON </files/5S_E_coli.json>` file.
```
:::

2. Now that you have an input file, navigate to the "Input/Output" tab to upload it.
3. Once the data is imported, you can continue work in several different ways:

* Download structure data to a different file type (Input/Output tab)
* Edit structure data (Edit tab)
* Add/remove/edit base pairs (Format tab)
* Add/remove annotations (Annotate tab)
* These options are explained in the “Tabs” and "Constraints" sections of this guide.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/IQhuC63lTRc" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

## Tabs

### The Input/Output tab

Here, you can:

* Upload input files
* Download output files

Together, these processes enable conversion between file formats.

#### How to upload input files

1. Click the "Open" button (Ctrl + O).
2. Select a file with one of the following supported input-file extensions:
   * `.xrna` -the native input-file format of XRNA (an extension of .xml)
   * `.xml` -a data-file format native to the internet. Highly similar to .xrna
   * `.json` -a general-purpose data-file format native to the internet. Logically equivalent to .xrna
   * `.str` -an image-file format containing nucleotide sequences and base pairs
   * `.svg` -an image-file format. Can be converted to other image-file formats (e.g. `.jpg`) using external tools
3. Upload the file.

#### How to download output files

Within the "Save File" row:

1. Provide an output-file name.
2. Select one of the following supported output-file extensions:
   * `.xrna` -the native input-file format of XRNA (an extension of .xml)
   * `.json` -a general-purpose data-file format native to the internet. Logically equivalent to .xrna
   * `.csv` -a comma-separated-values file. Logically equivalent to .xrna. Intended for use with RiboVision
   * `.bpseq` -a simplified, linear representation of the nucleotide sequence of the scene.
   * `.tr` -an XML-like format for R2DT. Contains nucleotides with 2D coordinates.
   * `.svg` -an image-file format. Can be converted to other image-file formats (e.g. .jpg) using external tools
3. Click the "Save" button (Ctrl + S)

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/y4h2paLWHkw" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

### The Viewport Tab

Here, you can manipulate the viewport with precision.
To clarify, the viewport is the main, right-hand component of XRNA. Within it, interactive 2D RNA diagrams are displayed.

#### Translation

"Translation" refers to the 2D displacement of the entire viewport. It is synchronized with click-and-drag within the viewport.

* `x`: The horizontal displacement. Rightward displacement is positive; leftward is negative.
* `y`: The vertical displacement. Upward displacement is positive; downward is negative.

#### Scale

"Scale" refers to the "in-and-out zoom" of the entire viewport. It is synchronized with the mousewheel within the viewport.

* Scrolling up with the mouse wheel increases zoom; this means zooming in.
* Scrolling down does the opposite.

#### Reset

In one step, the "Reset" button (`Ctrl + 0`) allows the user to reset all properties of the viewport:

* translate
* scale

This places all elements of the scene within view of the user.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/MIs1GuBH-gU" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

### The Edit Tab

Here, you can edit nucleotide data which has been uploaded and is displayed within the viewport.

1. Navigate within the viewport to the portion of the scene you plan to edit.
2. Next, do one of the following:
   * Select a constraint. Then, left-click on a nucleotide. Drag it to change relevant nucleotide-position data, according to the behavior of the constraint.
   * Middle-mouse-click and drag to select one or more nucleotide(s) or label(s). This will populate a menu which will allow you to precisely edit data, according to the behavior of an automatically determined constraint.
   * Select a constraint. Then, right-click on a nucleotide. This will populate a menu which will allow you to precisely edit data, according to the behavior of the constraint.

#### Notes

* Tutorials on how to use specific constraints are located in the "Constraints" section of this guide.
* Most constraints require you to click on a nucleotide that either is or is not base-paired to another nucleotide.
* Other constraints have more restrictive usage requirements.
* Annotations (i.e. labels) may be edited in the same exact way as nucleotides / nucleotide regions (described above).
* Annotations are not present in some input files, but they can be added within the "Annotate" tab.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/tHl9gmLTn7M" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

### The Format tab

Here, you can add, delete, or edit base pairs.

1. Navigate within the viewport to the portion of the scene you plan to format.
2. Do one of the following:
   * Middle-mouse-click and drag to select one or more nucleotide(s). This will populate a menu which will allow you to precisely format data, according to the behavior of an automatically determined constraint.
   * Select a constraint. Then, right-click on a nucleotide. This will populate a menu which will allow you to format base-pair data, according to the behavior of the constraint

#### Notes

* Tutorials on how to use specific constraints are located in the "Constraints" section of this guide.
* Most constraints require you to click on a nucleotide that either is or is not base-paired to another nucleotide.
* Other constraints have more restrictive usage requirements.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/2lVp4qd5kz8" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

### The Annotate Tab

Here, you can add, delete, or edit annotations (label lines or label content).

1. Navigate within the viewport to the portion of the scene you plan to annotate.
2. Do one of the following:
   * Select a constraint. Then, right-click on a nucleotide. This will populate a menu which will allow you to annotate nucleotides, according to the behavior of the constraint
   * Middle-mouse-click and drag to select one or more nucleotide(s). This will populate a menu which will allow you to annotate nucleotides, according to the behavior an automatically determined constraint

#### Notes

* Most constraints require you to click on a nucleotide that either is or is not base-paired to another nucleotide.
* Other constraints have more restrictive usage requirements.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/4C5YhlsAbWU" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

### The Setting Tab

Here, you can change settings which regulate how XRNA.js behaves.
Support for saving your settings is implemented by the pair of upload/download buttons.
Store the "xrna\_settings" file somewhere you will remember for later use.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/K9woUfHEBmg" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

### Constraints

#### The “Single Nucleotide” constraint

This constraint allows the user to interact with a single nucleotide.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/9wnXxSEwpKA" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “Single Nucleotide” constraint

This constraint allows the user to interact with a single base pair between two nucleotides.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/WJHwCsgauvA" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “RNA single strand” constraint

This constraint allows the user to interact with a contiguous series of single (non-basepaired) nucleotides.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/RTZIk_E-7SE" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The ”Single Helix” constraint

This constraint allows the user to interact with a contiguous series of base-paired nucleotides.
"Helix" is defined as two series of nucleotides which are mutually base-paired without gaps.
This involves at most two RNA molecules.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/hf2VUYS9clM" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “RNA sub-domain” constraint

This constraint allows the user to interact with a contiguous series of nucleotides constrainted by a helix.
This constraint involves only one RNA molecule.
This includes all nucleotides with nucleotide indices:

* Greater than or equal to the least nucleotide index in the helix
* Less than or equal to the greatest nucleotide index in the helix

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/sQ3vISkz1uc" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “RNA stacked helix” constraint

This constraint allows the user to interact with a contiguous series of helices between RNA molecules.
These helices may be separated by single (non-basepaired) nucleotides, but their mutually-basepaired status must resume outside the single-stranded regions.
This constraint involves at most two RNA molecules.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/UQu-zwaiKY8" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “RNA cycle” constraint

This constraint allows the user to interact with a cycle of nucleotides.
An RNA cycle is calculated as follows:

* Treat the nucleotides as nodes in a graph
* Treat the following as edges in the graph:
  + Base pairs between nucleotides
  + Neighboring nucleotides within a given RNA molecule (with nucleotide indices i, j such that abs(i - j) = 1)
* Calculate the smallest cycle within the graph involving the clicked-on nucleotide

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/2NSCglba3gg" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “RNA molecule” constraint

This constraint allows the user to interact with a single RNA molecule.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/oS2DwNiRZXg" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “RNA complex” constraint

This constraint allows the user to interact with a single RNA complex (which may contain more than one RNA molecule).

Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/AhODO78ol2A" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “Single Color” constraint

This constraint allows the user to interact with all nucleotides or labels with the same color

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/AhODO78ol2A" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### The “Entire Scene” constraint

This constraint allows the user to interact with all nucleotides or labels within the scene (i.e. viewport)

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/viejdl_aj8A" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

### Tips and Tricks

#### Hotkey combinations:

* `Ctrl + O` - opens a new input file
* `Ctrl + S` - saves a new output file, using the current output-file name and output-file extension.
  Fails if either of these are not provided.
* `Ctrl + 0` - resets the properties of the "Viewport" tab.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/evlD9BDW7U0" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

#### Middle Mouse click

In order to have a more enjoyable and intuitive experience with XRNA.js, try middle-mouse click-and-drag to highlight nucleotides.
This will:

1. Automatically change the constraint to an appropriate value which matches the nucleotides you selected
2. Populate a menu for precisely editing, formatting, or annotating data within the viewport.

#### Demo

```{eval-rst}
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; margin-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/pR54iebw1LI" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
```

## Contact Us

If you have any questions or feedback about XRNA-REACT, please get in touch with Anton S. Petrov at `anton.petrov@biology.gatech.edu`.
