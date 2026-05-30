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

I deployed openvscode-server in the Linux subsystem of the Chromebook. For convenience, I wrote a systemd service for openvscode-server and used the `autostart` plugin (`chromeos-autostart`) to make the Linux subsystem start automatically upon login.

The result: every time after boot, the host browser could not connect to `http://penguin:3000`. I had to manually enter the Linux subsystem and run `systemctl restart openvscode-server` before it became accessible.

### Root Cause

This issue is not actually related to openvscode-server itself. The Chromebook’s Linux subsystem (Terminal) follows a fixed startup sequence:

1. Start the container (penguin)
2. Load the user session
3. Execute the autostart plugin

However, when the autostart plugin triggers, systemd may not be fully ready yet, or network‑related services may not have finished initializing. As a result, although the service is marked as enabled, it does not actually start successfully (or the port is not yet bound to the container’s network interface).

### Solution

Add a restart command to the script executed by autostart. For example, I added the following line to my autostart script:

```bash
systemctl --user restart openvscode-server
```

Note that you must use `--user` because the openvscode-server service is a user‑level service, not a system‑level one. Without `--user`, systemctl will complain that the service cannot be found.

The complete autostart script looks roughly like this:

```bash
#!/bin/bash
# Start other services...
systemctl --user restart openvscode-server
```

After making this change and rebooting, the host can directly access `http://penguin:3000` without any issues.

## Bug 2: Black Screen and Color Blocks Caused by Snapdragon GPU Driver

### Problem Description

The HP Chromebook X11 uses the Snapdragon 7c Gen 2. The Adreno GPU driver on ChromeOS seems to have some minor issues, manifesting in two ways:

1. **Large black blocks in the terminal app**: Whether using the built‑in Terminal or crosh, scrolling often produces black squares/streaks that severely affect readability.

2. **Black screen when moving the mouse after opening opencode**: This is more serious. When I launch opencode (a TUI‑based AI programming tool) from the command line in the terminal, moving the mouse on the host desktop causes the screen to flash black. This is not occasional; it happens every time.

### Root Cause

Both problems are related to GPU‑accelerated 2D canvas rendering. ChromeOS enables GPU‑accelerated 2D rendering by default, but the Adreno GPU driver on the Snapdragon 7c Gen 2 may have compatibility issues with the 2D canvas implementation:

- The terminal app heavily uses canvas to render text and backgrounds. During scrolling, it constantly triggers redraws, and the driver cannot keep up, resulting in black blocks.
- opencode is a terminal TUI application whose rendering method triggers some kind of GPU state leak or deadlock, affecting the entire ChromeOS graphics stack. The redraw requests caused by mouse movement directly black out the display.

### Solution

Disable ChromeOS’s 2D canvas acceleration. In Chrome, navigate to:

```
chrome://flags/#disable-accelerated-2d-canvas
```

Set this flag to **Enabled** (note the name is `disable-accelerated-2d-canvas`; Enabled means acceleration is disabled), then restart the browser.

![disable-2d-canvas-flag](https://via.placeholder.com/600x200?text=disable-accelerated-2d-canvas)

After restarting:
- The black blocks in the terminal app disappear, and text renders normally.
- Moving the mouse after opening opencode no longer causes a black screen.

### Trade‑off

There is no free lunch. Disabling GPU‑accelerated 2D canvas means **all browser pages become noticeably slower**. Pages that heavily use canvas (such as Jupyter Notebook, VS Code Web, Figma, etc.) will feel laggy and have visible frame drops.

However, compared to a completely unusable black screen, a bit of lag is tolerable. I hope ChromeOS can fix this driver issue in the future, or that Qualcomm releases an updated GPU driver.

## Summary

Neither of these bugs is a major problem; once the cause is identified, the fix is straightforward. The first time you encounter them, though, they can be quite puzzling. The second black‑screen issue in particular took me several days to trace back to the 2D canvas acceleration.

The Chromebook X11 is a decent 2‑in‑1 device with good hardware, but it is still a niche platform. Users have to stumble upon and fill in these edge cases themselves. I’m documenting them here in the hope that it helps others facing the same issues.
