---
title: Two Bugs Encountered While Using Chromebook X11 and Their Solutions
date: 2026-05-30 20:00:00
tags:
  - Chromebook
  - ChromeOS
  - Bug
  - Snapdragon
  - GPU
categories: Others
---

Following up on the previous post {% post_link 从Fydetab-Duo到HP-Chromebook-X11-一个设备折腾党的自我修养 [From Fydetab Duo to HP Chromebook X11] %}, I have been using the HP Chromebook X11 for a while. The overall experience is quite good, but I encountered two annoying bugs. Here I record the solutions.

<!-- more -->

## Bug 1: Host Cannot Access Auto‑started Service in Linux Subsystem

### Problem Description

I deployed openvscode-server in the Linux subsystem of the Chromebook. For convenience, I wrote a systemd service for openvscode-server and used the `chromeos-autostart` plugin to make the Linux subsystem start automatically upon login.

The result: every time after boot, the host browser could not connect to `http://127.0.0.1:3000`. I had to manually enter the Linux subsystem and run `systemctl restart openvscode-server` before it became accessible.

### Root Cause

The root cause is completely unclear... I haven't seen any relevant discussion, and AI didn't find any similar issues either. Plus I never had this problem when using FydeOS.

### Solution

Add a restart command to the script executed by autostart. For example, I added the following line to my autostart script:

```bash
sudo systemctl restart openvscode-server
```

After making this change and rebooting, the host can directly access `http://127.0.0.1:3000` without any issues.

## Bug 2: Black Screen and Color Blocks Caused by Snapdragon GPU Driver

### Problem Description

The HP Chromebook X11 uses the Snapdragon 7c Gen 2. The Adreno GPU driver on ChromeOS for this SoC seems to have quite a few problems, manifesting in two ways:

1. **Large black blocks in the terminal app**: Whether using the built‑in Terminal or crosh, scrolling often produces black squares/streaks that severely affect readability.

2. **Black screen when moving the mouse after opening opencode**: This is more serious. When I launch opencode from the command line in the terminal, moving the mouse on the host desktop causes the screen to flash black. This is not occasional; it happens very easily.

### Root Cause

According to Gemini's answer, both problems are related to GPU‑accelerated 2D canvas rendering. ChromeOS enables GPU‑accelerated 2D rendering by default, but the Adreno GPU driver on the Snapdragon 7c Gen 2 may have compatibility issues with the 2D canvas implementation:

- The terminal app heavily uses canvas to render text and backgrounds. During scrolling, it constantly triggers redraws, and the driver cannot keep up, resulting in black blocks.
- opencode is a terminal TUI application whose rendering method triggers some kind of GPU state leak or deadlock, affecting the entire ChromeOS graphics stack. The redraw requests caused by mouse movement directly black out the display.

### Solution

Disable ChromeOS's 2D canvas acceleration. In Chrome, navigate to:

```
chrome://flags/#disable-accelerated-2d-canvas
```

Set this flag to **Disable**

After restarting:
- The black blocks in the terminal app disappear, and text renders normally.
- Moving the mouse after opening opencode no longer causes a black screen.

### Trade‑off

This feature was originally for GPU-accelerated 2D graphics rendering. After disabling it, the related rendering will inevitably become slower...

But to be honest, compared to a completely unusable black screen, a bit of lag is at least tolerable. I hope ChromeOS can fix this driver issue in the future (although it's not very realistic—it seems 7c devices have always had various problems).

## Postscript

Choosing commercial products was originally for better usability and fewer problems... Even Chromebooks still have issues... The constant small problems with Fydetab Duo before suddenly seem more acceptable...
