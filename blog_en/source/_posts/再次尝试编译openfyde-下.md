---
title: Trying to compile openfyde again (Part 2)
categories: Others
date: 2025-06-19 23:09:23
tags: ['fydeos', 'openfyde']
---

After changing the motherboard and CPU, everything went smoothly... I'm not sure whether it was the motherboard or CPU that was causing the issue before - compilation would fail after tens of minutes, and subsequent attempts would fail at different points. This was actually the first time in my life I encountered real hardware incompatibility issues... Another "first time in my life" experience.

So here comes the second part of compiling openfyde. Unfortunately, during the time I was struggling with the faulty hardware, the prebuilt image for r132-dev version was already updated... Missed the hot topic...

<!-- more -->

Actually, I had already reached the final steps last time, but got stuck when compiling the massive Chrome project due to my hardware issues. After completing the compilation, we can generate the image that can be directly flashed to a USB drive.

This time I actually solved a problem from last time - how to get a flashable image using Rockchip tools after obtaining `chromiumos.bin`. The tutorial gives us a bin image, while official prebuilt versions are in img format, clearly different. Also, unpacking official openfyde prebuilt images with Rockchip tools seems to yield some content not included in the bin.

Actually, the official tools are available - it's [openFyde/rk3588-image-maker](https://github.com/openFyde/rk3588-image-maker). But last time as a complete beginner, I couldn't understand what was happening during compilation, let alone figure out what all the projects under openfyde were for.

The steps to generate the required image are simple: first clone this project, then go to the root directory and run the following commands as instructed:

```bash
# Mount required content from bin
./map_chromiumos_image.sh /PATH/TO/CHROMIUMOS/IMAGE.bin --board fydetab_duo
# Generate update.img
./rk3588-mkupdate.sh
```

Note: For the r132 version I compiled this time, the `Image/parameter.txt` file in the project needs modification - change `0xa0006d@0x00b02040(STATE)` to `0xa0006d@0x00b02040(STATE)`.

Then copy the generated image and flash it using the official Rockchip tool. The compiled image this time has working touchscreen, Bluetooth, and Linux containers - should be problem-free.

Next step is to try adding stuff to it!
