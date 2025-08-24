---
title: Flashing Mokee on My Nutro Pro
categories: Others
date: 2021-02-12 20:56:03
tags: ['phone', 'Mokee', 'android']
---

I've had my Nutro Pro for over three years now, and it's been slowing down significantly. The storage space is also running low, making it hard to use comfortably. So, I decided to give it one last shot by flashing a new ROM. If it works out well, I'll continue using it for another year or so. Otherwise, I'll have no choice but to buy a new phone with my end-of-year bonus.

## Tools Needed

Speaking of the last time I flashed my phone, it was quite some time ago. Back then, custom ROMs were all the rage as people not only enjoyed their phones but also played around with different operating systems. Nowadays, unless you have a specific need or your original system is really slow and unresponsive, most people probably won't consider flashing their phones. Flashing has become more complicated than it used to be; in fact, I needed a special cable for this task. Fortunately, there are people who compile and maintain ROMs for this device. Here's what you'll need:

1. A 9008 flash cable (about 10 RMB on Taobao, with a button that stays pressed when clicked and releases when clicked again)
2. [TWRP for Nutro Pro](https://download.mokeedev.com/odin.html)
3. Flashing tools (`QPST` for Windows, `qdl` for Linux; you can install `qdl` using `yay -S qdl`, which is an AUR package)
4. The base and system packages for Nutro Pro (you can find them at [MokeeDev](https://download.mokeedev.com/odin.html))

## Flashing Steps

The flashing process follows the instructions provided by the developer, but since I'm not a Windows user, I used different tools than what was recommended. Additionally, due to inconsistencies between the cable I bought and the one mentioned in the tutorial, the actual method of entering 9008 mode was slightly different (the button on my cable stays pressed for a moment before releasing).

As for backing up files before flashing... I only backed up my WeChat chat records using the official WeChat backup feature. Other data can be left as is since I don't use phone calls or SMS messages much anymore. My phone was already running the latest official system at the time of flashing, and I hadn't performed any previous flashes.

### Place Base and System Packages on SD Card
I have a 32GB SD card in my phone, so I placed both the base and system packages directly on it. It should be possible to place them on the internal storage as well since there's an option during flashing, but I haven't tried it yet.

### Run qdl on Computer and Wait for Phone Insertion
Download TWRP and install `qdl`. Extract the TWRP files and navigate to the extracted directory. Run the following command to see "Waiting for QDL tty...":

```
sudo qdl prog_emmc_firehose_8953_ddr.mbn rawprogram_unsparse.xml
```

### Enter 9008 Mode on Phone and Start Flashing

1. Long-press the power button to turn off the phone.
2. Connect the 9008 cable to the USB port of your computer.
3. Simultaneously press the volume up and down keys.
4. Keep the 9008 cable button unpressed while inserting it into the phone. The top-left indicator light on the front will briefly flash red; once you see this, immediately press the button. After that, the indicator light will keep flashing red, indicating that the phone is in 9008 mode. You can release your fingers from the volume keys at this point.
5. If everything goes smoothly, you should see the TWRP being flashed on your computer terminal. Wait for it to finish.

### Flashing Mokee

Follow the instructions provided by [MokeeDev](https://bbs.mokeedev.com/t/topic/14503) as closely as possible. It took me about 5-10 minutes to complete the flashing process. After rebooting into the system and making necessary settings, you can start using it.

## Short-term Usage Experience/Minor Issues

I've been using it for two days so far, and I feel that it's significantly more responsive than the original system, though not as good as expected. The system-level performance is excellent, but app performance still has some issues. It seems that the CPU can't keep up with the performance demands of domestic apps. However, since apps can be moved to the SD card, it should last me another year without issues. The current problems I've encountered are:

1. Deleting fingerprints requires renaming the fingerprint first before deletion.
2. Force stopping processes is very aggressive; if you don't interact with the phone for a while, non-background processes will crash. Some comments suggest that this is due to domestic vendors not following Google's requirements when developing apps, which is a common issue for apps similar to the original Android system.
3. Occasional small lags, especially when Bilibili has many弹幕.

## Afterword

The system-level performance after flashing is excellent, but there are still occasional issues within commonly used apps. Then I remembered the "evolution" process of WeChat...

So, I tried replacing it with an older version... It worked! Bilibili no longer lags, and the startup speed of Alipay is noticeably faster.

After realizing that newer versions of apps might be causing performance issues, I tried replacing them with older versions... It worked! Bilibili no longer lags, and the startup speed of Alipay is noticeably faster. This makes me think about how apps have evolved over time, often becoming more resource-intensive.

Well... what can I say... Here's a list of historical versions of the app that might help others facing similar performance issues...
