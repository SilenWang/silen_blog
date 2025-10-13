---
title: TypeScript和JavaScript中的调试技巧
categories: Coding
date: 2025-10-13 23:08:14
tags: ['调试', 'Typescript', 'Javascript', '前端开发']
---

正如我同学所说，我现在有一点"被迫转码"的趋势... 目前维护的几个网站都有各自的前端、后端数据库，部分还有测试代码和迁移代码，涉及语言从生物信息学常用的Python到之前几乎不用的JS、TS、HTML、C#等。实际工作中，不可能有太多时间从头系统学习每门语言，因此掌握最基本的调试方法至关重要。其中JavaScript和TypeScript的调试方式，我觉得特别值得记录，因为它们与其他语言相比确实有些独特之处。

<!-- more -->

## Web前端框架中的调试环境

与其说是纯JS/TS调试，不如说是Web前端框架环境下的调试更为准确。现代前端框架（如Vue、React、Angular）与传统的HTML文件结构既有相似之处，也有明显差异。这些框架本质上是通过各种编译和构建工具，将模板、脚本和样式代码转换为最终的网页代码，再由浏览器渲染显示。

这种混合编程的模式类似于Snakemake的工作流文件，单个文件中可能包含多种语言元素。在实际开发中，我们主要需要调试的是模板部分（HTML-like）和脚本部分（JavaScript/TypeScript）。

## 模板部分的调试技巧

### 1. 直接显示变量值
在Vue模板中，可以使用双花括号语法直接显示变量：
```vue
<p>{{ value }}</p>
```

在React中，则使用单花括号：
```jsx
<p>{value}</p>
```

### 2. 条件性调试显示
有时我们只想在开发阶段显示调试信息：
```vue
<template>
  <div>
    <p v-if="isDevelopment">调试信息: {{ debugValue }}</p>
    <!-- 正常内容 -->
  </div>
</template>
```

## 脚本部分的调试方法

### 1. 基础console调试
```typescript
console.log('变量值:', variable);
console.table(arrayData); // 以表格形式显示数组或对象
console.dir(object); // 显示对象的属性
```

### 2. 分组日志
```javascript
console.group('用户信息');
console.log('姓名:', user.name);
console.log('邮箱:', user.email);
console.groupEnd();
```

### 3. 条件断点和debugger
```typescript
if (someCondition) {
  debugger; // 浏览器会在该行自动暂停
}
```