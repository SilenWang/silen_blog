---
title: Flashing Mokee ROM on My Nuts Pro
categories: Daily
date: 2021-02-12 20:56:03
tags: ['phone', 'Mokee', 'android']
---

I've been using my Nuts Pro for over three years, and everything is slow, and the storage is always running out. It's time to make a final effort... Desperate times call for desperate measures, so I'm going to try flashing it. If it works well, I can keep using it for another year and a half. If not, I'll have to buy a new one with my end-of-year bonus.
<!-- Summary -->
<!-- more -->

## Preparation

Talking about the last time I flashed my phone, that was a long time ago... Back then, Android phones were just starting to become popular, and people liked playing with different custom ROMs. Now, unless you have a special need or can't stand your original system anymore, most people probably don't flash their phones anymore. Flashing has also become more complicated now... At least my phone needs a specific cable. Of course, I'm grateful to those who are compiling and maintaining the ROM for this device. So, here's what you'll need:

1. A 9008 flashing cable (available on Taobao for 10 yuan)
2. [TWRP for Nuts Pro]()
3. Tool for flashing TWRP (`QPST` for Windows, `qdl` for Linux, install with `yay -S qdl`, it's an AUR package)
4. The corresponding [base package](https://download.mokeedev.com/odin.html) and [system package](https://download.mokeedev.com/odin.html) for Nuts Pro (`RADIO-odin*` and `MK90.0-odin*`)

## Steps

The steps are basically following the developer's [tutorial](https://bbs.mokeedev.com/t/topic/14503), but since I'm not a Windows user, the flashing tool is different (although I tried both, and the Windows tool didn't work for me, it only confirmed that my phone was in 9008 mode). Also, because of the inconsistency between the cable I bought and what's mentioned in the tutorial, the actual way to enter 9008 mode is not consistent with the tutorial (the button on my cable stays pressed when you press it once, and then pops back when you press it again).

As for backing up files before flashing... I only backed up my WeChat chat records. The rest can be ignored... After all, I haven't been making phone calls or sending text messages in a long time... Also, my phone was updated to the latest official system before flashing, and there had never been any previous flashing.

### Put base package and system package on SD card
My phone has a 32GB SD card. I put both the base package and the system package directly on the SD card. It should be possible to put them in internal storage as well, since I saw an option during the flashing process. However, since I haven't tried it yet, I'm not sure.

### Run qdl on your computer and wait for the phone to connect
Download TWRP, install `qdl`, unzip it, enter the unzipped directory, and run the following command until you see `Waiting for QDL tty...`. If it fails, try running it again.

```
sudo qdl prog_emmc_firehose_8953_ddr.mbn rawprogram_unsparse.xml
```

### Put phone in 9008 mode and start flashing

1. Long press the power button to shut down.
2. Connect the 9008 cable to your computer's USB port.
3. Press the volume + and - buttons simultaneously.
4. Keep the 9008 cable button in an unpressed state, then insert it into the phone. The red light on the front left corner will briefly flash. When you see the red light, press the button immediately. After that, the phone's light will keep flashing red, indicating that it has entered 9008 mode. You can release the volume key at this point.
5. If everything goes smoothly, you should see TWRP being flashed on your computer terminal. Wait for it to finish.

### Flash Mokee ROM

Follow the [tutorial](https://bbs.mokeedev.com/t/topic/14503) completely. I did exactly what was written without doing anything extra... It took about 5-10 minutes. After that, restart and set up the system as needed.

## Short-term experience/problems

After using it for two days, I feel that it's much smoother than the original system most of the time, but not as good as expected... It seems that the CPU can't keep up with the performance demands of domestic apps. However, since apps can be moved to the SD card, it should be fine for another year. The problems I've encountered so far:

1. When deleting fingerprints, you need to rename them first before deleting.
2. Processes are killed very aggressively. If you don't operate for a while, non-background processes will definitely hang. Of course, this is what many comments say... It's the result of domestic vendors not following Google's requirements when developing apps... A common problem with Android-like systems using domestic apps.
3. Small lags still exist, especially when there are many Bilibili danmu.

## Postscript

After flashing, the system level is very smooth, but occasionally there are problems within commonly used apps... Then I suddenly thought about WeChat's "evolution"... Hmm...

So I tried replacing old versions of the app... Well, that would be great. Watching videos on Bilibili no longer lags, and the startup speed of Alipay has visibly increased.

What can I say... I'll write a list of historical app versions later...
