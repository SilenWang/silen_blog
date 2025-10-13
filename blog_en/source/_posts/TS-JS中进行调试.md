---
title: Debugging Techniques in TypeScript and JavaScript
categories: Coding
date: 2025-10-13 23:08:14
tags: ['Debugging', 'Typescript', 'Javascript', 'Frontend Development']
---

As my classmate said, I'm currently experiencing a situation of "being forced to code"... The several websites I maintain each have their own frontend, backend databases, and some even have test code and migration code, involving languages ranging from Python commonly used in bioinformatics to JS, TS, HTML, C#, which I rarely used before. In actual work, there isn't much time to systematically learn each language from scratch, so mastering the most basic debugging methods is crucial. Among them, I find the debugging methods for JavaScript and TypeScript particularly worth recording, as they do have some unique aspects compared to other languages.

<!-- more -->

## Debugging Environment in Web Frontend Frameworks

It's more accurate to say it's debugging in a web frontend framework environment rather than pure JS/TS debugging. Modern frontend frameworks (such as Vue, React, Angular) have both similarities and obvious differences compared to traditional HTML file structures. Essentially, these frameworks use various compilation and build tools to convert template, script, and style code into final webpage code, which is then rendered and displayed by the browser.

This mixed programming model is similar to Snakemake workflow files, where a single file may contain multiple language elements. In actual development, what we mainly need to debug are the template part (HTML-like) and the script part (JavaScript/TypeScript).

## Debugging Techniques for Template Parts

### 1. Directly Display Variable Values
In Vue templates, you can use double curly brace syntax to directly display variables:
```vue
<p>{{ value }}</p>
```

In React, use single curly braces:
```jsx
<p>{value}</p>
```

### 2. Conditional Debug Display
Sometimes we only want to display debug information during the development phase:
```vue
<template>
  <div>
    <p v-if="isDevelopment">Debug Info: {{ debugValue }}</p>
    <!-- Normal content -->
  </div>
</template>
```

## Debugging Methods for Script Parts

### 1. Basic Console Debugging
```typescript
console.log('Variable value:', variable);
console.table(arrayData); // Display arrays or objects in table format
console.dir(object); // Display object properties
```

### 2. Grouped Logs
```javascript
console.group('User Information');
console.log('Name:', user.name);
console.log('Email:', user.email);
console.groupEnd();
```

### 3. Conditional Breakpoints and Debugger
```typescript
if (someCondition) {
  debugger; // Browser will automatically pause at this line
}
