R2DT Compare — quick start for curators
=======================================

This folder starts a private tool on your own computer so you can compare
RNA structures, edit base pairs, and save your work locally.

You do not need to install R2DT, Python, or Git.


Before you begin
----------------
Install Docker Desktop (free) and open it once so it is running:

  https://www.docker.com/products/docker-desktop/

Wait until Docker says it is running (whale icon idle / no “starting…” message).

Linux users: install Docker for your distribution instead, then make sure you
can use Docker without special admin steps each time. Ask IT if unsure.


Start
-----
Pick the file for your computer:

  Mac     →  Start-macOS.command   (see “Mac security warning” below — do this first)
  Windows →  Start-Windows.bat     (double-click)
  Linux   →  Start-Linux.sh
             (if a double-click does nothing, open a Terminal in this folder
              and type:  ./Start-Linux.sh )

A Terminal / Command window will open. Leave it open while you work.


Mac security warning (expected)
-------------------------------
macOS often blocks the starter the first time with a message like:

  “Apple could not verify Start-macOS.command is free of malware…”

That is normal for an unsigned script you downloaded. It is not a virus scan
of R2DT itself. Clear it once like this:

  1. In Finder, go to this folder.
  2. Control-click (or right-click) Start-macOS.command — do not double-click yet.
  3. Choose Open.
  4. In the dialog, click Open again.

After that, double-click works as usual.

If macOS only offers “Done” / “Move to Bin”, open System Settings →
Privacy & Security, scroll to the message about Start-macOS.command, and
click Open Anyway. Then try Open again from Finder.

Alternative (Terminal), from this folder:

  xattr -d com.apple.quarantine Start-macOS.command
  ./Start-macOS.command


First time
----------
The first start downloads the R2DT toolkit. That can take several minutes and
needs a network connection. Later starts are much faster.


What you should see
-------------------
Your browser should open the Compare page automatically:

  http://127.0.0.1:8765/compare

If it does not, copy that address into Chrome, Firefox, Safari, or Edge.

From there you can upload a reference and a model, choose matching chains,
run a comparison, edit base pairs, and export results.


Your work stays on this computer
--------------------------------
Jobs and edits are stored in a folder named .r2dt-workstation in your home
directory. Closing the browser does not delete them. Updating or restarting
keeps that folder.


When you are finished
---------------------
Click the Terminal / Command window and press Ctrl+C to stop the tool.
You can start again any time by double-clicking the same Start file.


If something goes wrong
-----------------------
• “Docker is not running” or “Docker was not found”
    Open Docker Desktop, wait until it is ready, then try Start again.

• Mac: “Apple could not verify … malware” / cannot be opened
    Control-click Start-macOS.command → Open → Open.
    Or: System Settings → Privacy & Security → Open Anyway.
    Or in Terminal (in this folder):
      xattr -d com.apple.quarantine Start-macOS.command
      ./Start-macOS.command

• Linux says “Permission denied”
    In Terminal, in this folder, run:  chmod +x Start-Linux.sh
    then try ./Start-Linux.sh again.

• Port already in use / page already open
    The tool may already be running — just open the Compare link above.

Need help? Contact the person who sent you this folder and mention any
exact message shown in the Terminal window.

(This pack is also linked from the Local workstation docs page as
r2dt-workstation-start.zip.)
