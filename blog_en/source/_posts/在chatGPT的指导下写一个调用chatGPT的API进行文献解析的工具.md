---
title: In the Guidance of chatGPT, Write a Tool to Parse Literature Using chatGPT's API
tags:
  - chatGPT
  - AI
categories: Coding
date: 2023-03-27 01:26:48
---


This title is a bit too nested... but the fact is, I suddenly felt capable! I finally have a chance to become a full-stack developer! To showcase my capabilities enhanced by chatGPT, I decided to develop a tool that uses chatGPT to read papers (I've also considered using it for meta-analysis).

<!-- Abstract part -->
<!-- more -->

## Start Digging

About two weeks ago, I had the idea of doing this. Previously, I registered an account and tried using it to assist in writing a static service with token verification. Later, I saw my colleague share [chatPDF](https://www.chatpdf.com) and mentioned that he uses chatGPT to quickly read papers and write reviews. At that time, I recalled the pain of reading papers here over three years and thought about developing a tool to accelerate my reading using chatGPT to assist in my work, which could also serve as a demo for popularity and potentially gain several hundred stars, bringing me closer to becoming an independent full-stack developer who has always wanted. Plus, before, I couldn't do it because learning everything from scratch was too time-consuming and lacked the ability. With chatGPT's assistance, it should be much easier.

## Various Pitfalls

Then, let's get down to business! (Three-minute热度 is most擅长 this), I thought that with such an external aid, developing a simple demo wouldn't take more than three days? So, over the past two weeks, I've been going through various pitfalls...

### Paying Money Is So Difficult...

I didn't expect the first pit to be paying money... After registering for a new OpenAI account, I received 30 dollars of free API credits for trial use. However, it seems that too many users have registered, and they are reducing or even stopping the free credits... When I registered, my account had no credits at all, so I couldn't use the core API functionalities. I had to pay out of pocket... But it's not as simple as that; OpenAI does not provide services for mainland China developers... Therefore, my Visa card couldn't be registered... I had to apply for a virtual credit card... Then this virtual credit card... I'm not sure if it's a legitimate virtual currency... The largest Depay requires using a certain type of virtual currency to recharge??? Virtual currencies make me feel scared; I dare not mess with them. In the end, I found a large amount that I could pay directly (I might buy塞尔达 later), and finally solved the problem (it seems it was still resolved by paying money). After topping up, I started developing!

emmmmmm... Then after half a day...

### Limitations of chatGPT

chatGPT can execute various tasks, but this doesn't change the fact that it is a text generation AI. Its fundamental function is to generate fluent and error-free natural language text. Initially, when using it, I didn't encounter many problems (perhaps because the questions were simple), but as usage frequency increased, the frequency and severity of encountering issues also rose significantly. In general, I feel that due to asking overly broad questions or asking questions that may not have answers, I receive various irrelevant answers from chatGPT. This leads to two issues:

1. During actual development, it wasn't as fast as expected because I adopted a very new development framework (`Pynecone`). Its introduction is beautiful; once developed, frontend/backend/API are all included, and multi-platform quick deployment is possible. However, due to its immaturity, chatGPT couldn't give me any suggestions. I had to learn the documentation and official examples to complete necessary functional implementations. Moreover, it took a lot of time on some very basic functionalities...

2. Making chatGPT execute specific tasks smoothly isn't easy... And I'm using `gpt-3.5-turbo` (which didn't have version 4 at the start; now I can't queue up), so it's not as easy to specify the return answer format with gpt4 (I want to use JSON format for easier parsing). Therefore, "programming towards prompts" is no joke... To use it well, debugging prompts is still very necessary.

In summary, chatGPT is strong, but... even gpt4, I believe, isn't as powerful as the Microsoft promotional video suggests... To use it smoothly, you need to have some usage skills and become shaped like chatGPT (misunderstanding).

### Switching Development Frameworks

Continuing from the previous text... `Pynecone` really isn't mature enough... After a week and a half, I couldn't take it anymore and switched to `Gradio`... How to say it? The tool and framework's Matthew effect shouldn't be violated casually... Especially since my intention was to gain popularity, finding the fastest implementation method is more important than overthinking.

### SinGAN you sin!

Although just a joke, I really didn't expect... I got stuck for an entire afternoon and one night while deploying the demo on HuggingFace. HuggingFace is a great platform (a space provides 8 cores and 16GB, what a charitable oil tycoon), but not allowing a new project to be completely empty really made it difficult for me as a beginner...

I experienced manual conflict resolution, forced code pushes, and ultimately compromised by separating the original project's code from HuggingFace. I really couldn't get any effective help on this issue from chatGPT.

### Not Enough Time

Worse still, there was rain on my parade... (The road to becoming a national treasure is irreversible)

## Help Provided by chatGPT

Although I've complained that chatGPT isn't as good as it seems, the benefits it brings are real. In this project, the help provided by chatGPT includes but is not limited to:

1. Helping me understand similar open-source code quickly, allowing me to quickly grasp other implementation methods.
2. Answering various usage questions for Gradio (a mature tool can provide many practical suggestions).
3. Assisting in generating English README (although no one reads it).
4. Increasing my confidence that "I can do it!"

## Feeling the Pain of Independent Development

These two weeks also deeply made me feel the pain of independent developers or personal自媒体. Similarly, as a literature reading assistant, the [chatPaper](https://github.com/kaixindelele/ChatPaper) project and my content both started disclosing code on GitHub two weeks ago (my case was uploading code incrementally...), but now it has over 4000 stars... I...

This feeling is probably similar to the anxiety of starting a自媒体, only to find that no one watches your videos.

## Project Address

However, regardless of this, the purpose of my project is to develop something useful that I personally feel good about. If popularity can't be gained... then it's fine... My project is called [ResearchGPT](https://github.com/SilenWang/ReviewGPT), and the demo is hosted on Huggingface at [Demo](https://huggingface.co/spaces/SilenWang/ReviewGPT). Those interested can try it out~ Of course, this blog post also embedded the demo.

<div style="display:flex;justify-content:center;align-items:center;overflow:scroll">
   <iframe
      	src="https://silenwang-reviewgpt.hf.space"
      	frameborder="0"
      	width="850"
      	height="700"
   ></iframe>
</div>
