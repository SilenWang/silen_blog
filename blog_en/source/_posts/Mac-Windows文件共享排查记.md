---
title: Troubleshooting File Sharing Between Mac and Windows — and a Rant About Windows' Messy Design
categories: Others
date: 2026-08-10 12:00:00
tags: ['macOS', 'Windows', 'SMB', 'File Sharing', 'Network Troubleshooting', 'Win11', 'rant']
---

I wanted to share files on my company's Mac Mini with other Windows computers over SMB. The Mac clearly had "File Sharing" turned on, yet the Windows machines refused to see the machine, let alone the shared folders. After half a day of tinkering, I finally found the real culprit — with the help of a pile of command-line tools. This post records the whole troubleshooting process, and also gives a good rant about Windows' junk-pile design that has only gotten messier from Win10 to Win11.

<!-- more -->

Let me state the conclusion first: nine times out of ten, this kind of problem comes down to things not being aligned on the two sides — either Windows "Network Discovery" is off, or the SMB protocol versions don't match, or the credentials aren't filled in correctly. But the truly infuriating part is that these settings are scattered across several places on Windows, and just hunting them down is enough to drive you crazy.

## First, confirm the network is actually up: direct IP connection

The first step isn't to fiddle with those flashy toggles, but to confirm the underlying network actually works. It's simple: run `ipconfig` on Windows to get the IPv4 address, then on the Mac press `Command+K` in Finder and enter `smb://IP`.

```
ipconfig        # Get the IPv4 address on Windows
# On Mac Finder: Cmd+K → smb://192.168.x.x
```

If it connects, the problem is basically in the "Network Discovery" layer; if it doesn't, you need to dig deeper into SMB protocol and credentials.

To access in the reverse direction, just type `\\Mac's-IP` in the Windows File Explorer address bar. Note the two leading backslashes, e.g. `\\192.168.1.5`. This bypasses Network Discovery and connects directly to the IP — if this works, it means the Mac's sharing service itself is fine.

## When Windows can't connect to `\\IP`

Of course, if even `\\IP` fails, you should work through the following checks one by one:

**First check on the Mac side:**

1. The Mac user must have a login password set — SMB usually doesn't allow accounts with empty passwords.
2. System Settings > General > Sharing > File Sharing > "Options...", and confirm the Mac username you plan to log in with is checked.
3. Confirm the shared folder's permission for that user (or Everyone) is "Read & Write" or "Read Only".
4. Use the username format `MacComputerName\username`, e.g. `MacBook-Pro\zhangsan`.

**Digging deeper on the Windows side:**

```powershell
ping MacIP
telnet MacIP 445                       # Or Test-NetConnection -ComputerName <MacIP> -Port 445
sc.exe query lanmanworkstation
net use \\MacIP /user:username password   # Returns specific error codes like 1312/1326
secpol.msc                              # Check whether "Network access: Sharing and security model for local accounts" is "Classic"
```

The error codes from `net use` are very useful — numbers like 1312/1326 directly tell you which layer the problem is in.

**Digging deeper on the Mac side:**

```bash
sudo lsof -i :445   # Confirm smbd is listening
# Restart the SMB service:
sudo launchctl unload -w /System/Library/LaunchDaemons/com.apple.smbd.plist
sudo launchctl load -w /System/Library/LaunchDaemons/com.apple.smbd.plist
# Loosen the protocol negotiation versions:
sudo defaults write /Library/Preferences/com.apple.smbd.plist SmbMinVersion "SMB_2_00"
sudo defaults write /Library/Preferences/com.apple.smbd.plist SmbMaxVersion "SMB_3_00"
```

The "Block all incoming connections" option in the Mac firewall should be turned off, and you can also try accessing via `\\ComputerName.local` instead of the IP.

Two key pieces of evidence came out at this point, which narrowed the problem down directly:

1. `sudo lsof -i :445` showed the `smbd` process listening on port 445, on both IPv4 and IPv6 — **the Mac's SMB service is completely normal**.
2. The `Test-NetConnection` result on Windows: `TcpTestSucceeded : True` (RemoteAddress 172.16.1.25, SourceAddress 172.16.1.21, InterfaceAlias WLAN) — **port 445 from Windows to the Mac is fully open, no problem at all at the network layer**.

With the network working and the service normal, the only remaining possibilities are SMB protocol version negotiation and authentication/credentials.

## The final solution

In the end, two commands run in an elevated Windows PowerShell did the trick:

```powershell
Set-SmbClientConfiguration -RequireSecuritySignature $false -Force
Set-SmbClientConfiguration -EnableInsecureGuestLogons $true -Force
```

Right after running them, it connected. These two commands disable the SMB signing requirement and allow insecure guest logons — meaning the real blocker was in the protocol negotiation layer, not the network or the firewall.

## The rant: Windows' junk-pile design

This whole ordeal gave me a real taste of how chaotic Windows' design is. Going from Win10 to Win11, it hasn't gotten more consolidated — it's gotten more fragmented:

**First, there are two shells.** cmd and PowerShell. The same setting sometimes has to be changed in cmd, sometimes in PowerShell, and the two have completely different syntaxes. Even a configuration command like `Set-SmbClientConfiguration` this time had to be run in an *elevated* PowerShell — the average user can't even get their foot in the door.

**Second, there are two settings UIs.** "Control Panel" and "Settings" coexist, and the same setting exists in both places but under different paths and names. Network sharing settings hop between the Control Panel's "Network and Sharing Center" and the new Settings app, and their location shuffles around with every version update. Tutorials from ten years ago point at Control Panel, this year's tutorials are all about Settings — you have no idea which one you should be looking at.

**Third, the same thing bounces around between multiple places.** Take the SMB signing toggle alone — it hops back and forth between Control Panel, Settings, and PowerShell commands. To troubleshoot one problem, you have to be able to guess UI paths, memorize command lines, and know your way around policy items at the same time — this is not a system designed for humans, it's clearly built for veteran devs who've already survived the handover.

**The macOS side, by contrast**, keeps everything neatly concentrated: System Settings > General > Sharing, one chain that handles it all. File sharing, options, and permissions are all on that single page — a couple of clicks and you're done. Put the two side by side, and the verdict writes itself.

## Summary

The trap this time can be summed up in one sentence: **the network was fine, the service was fine, the whole problem was in Windows' shoddily-assembled protocol settings**. If you also run into Windows not seeing your Mac, first rule out network issues with a direct IP connection, then check Network Discovery and the workgroup, and finally use commands like `Test-NetConnection` and `lsof -i :445` to lock down the scope before applying the fix. As for Windows' junk-pile design... you just have to adapt to it.
