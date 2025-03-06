---
title: Rewriting My Project in Go Language
categories: Script
date: 2019-06-02 14:46:49
tags: ["golang"]
---

To learn Go language, I rewrote a script that was used frequently in my previous work. Compared to Python, which is simple, quick, easy to understand and has many useful third-party modules with abundant learning materials in both Chinese and English, if it wasn't for performance reasons, I really wouldn't want to switch...
<!-- Abstract part -->
<!-- more -->

## A Go Language Beginner's View
Since most of my larger programs were written in Python and R, which are not strong-typed languages, the difficulty mainly lies in using third-party modules and handling logical bugs designed by myself.

Go is different... it is a strongly typed language. Variables must be declared before use, and once declared, their type cannot be changed later. This made me feel uncomfortable for a long time when I first started learning Go...

Although Go provides interfaces to some extent to increase usability, as a beginner who didn't have similar concepts in Python, R, or Bash, it took me quite some time to understand what an interface actually is.

After solving the basic variable usage issues, there were differences in error handling methods, fewer third-party modules, and lack of Magic writing techniques, among other things. These problems had to be solved one by one... ultimately leading to this rewrite taking two whole weeks to initially complete...

Although it was difficult, once you experience its speed, many people would be as excited as me (although during this rewrite, I also optimized the algorithm, but that optimization is absolutely not enough to bring a 200-fold difference).

Well, the world changes so fast, learning something new is never wrong~

Here's the address of the rewritten project: [FUEX](https://github.com/SilenWang/FUEX)

## Quickly Mastering Go with Python as Reference
Just like how I always found it hard to think in English when learning it from Chinese, writing this practice project was inevitably influenced by my way of thinking in Python. However, Python and Go are two very different languages, some things that are used in Python have similar or even no equivalents in Go, so here's a summary of the alternative implementation schemes I encountered during this development.

### Variables/Data Types

#### Variable Declaration and Initialization
In Python, there is no separate step for variable declaration and initialization. You can directly use `var=value` to use it, and the type of the variable can be changed later when processing. For example:

```python
var1 = True
var2 = 1
var1 = var2 if var2 else True
```

In the above code, `var1` and `var2` are automatically assigned types (boolean, int), and although their types are different, you can pass the value of `var2` to `var1`, and after passing, `var1` will change from boolean to int.

While all variables need to be declared like this, Go provides a simple way: `var1 := 1`. This allows you to declare and assign a variable at the same time, with the type of the variable automatically specified. However, since it's automatically specified, if there are special requirements for the type, it might still be better to declare it manually.

#### Go is a Strongly Typed Language

The term "strongly typed" refers to the detailed classification of variables in the language and the strict restrictions on passing and using variables of different types during program execution (personal understanding). In Python, the variable types are numbers, strings, logic, and empty (`None`), with complex ones being numbers further divided into floating-point and integers. Even Python has gotten rid of constants, only having variables, but it recommends using all uppercase letters to indicate constants in language writing suggestions (actually, the value can still be changed).

In Go, there are constants (`Const`) and variables (`Var`). While variables have integer, float, logic types among others, integers and floats are further divided into different lengths such as `int16`, `int32`, `int64`. Even the `uint` I'm not quite clear about... These variable types cannot be automatically converted during use. For example, if a function gives you an `int16`, and the next function needs an `int64`, you still need to manually convert it once... These type conversions are really very uncomfortable for me who is used to Python.

Another thing to note is that in Go, `byte` and `string` are two different types. Characters use single quotes ('), while strings use double quotes (").

#### Complex Data Types----Arrays/Slices/Collections/Structures
Just like Python has Tuple, List, Dict, Set for storing multiple elements or mapping data, Go also has data types that meet the needs of data storage and processing:

Array: Arrays are similar to Lists or Tuples. First, the length of an array in Go is fixed and cannot be changed once declared. At the same time, although the contents stored in an array can be mutable, the type of the stored elements must be fixed.

Slice: Similar to arrays, but its length is variable.

Map: Equivalent to Dict, similarly, both the Key and Value types are fixed and unchangeable. The Keys used by Map can be any type that supports logical comparison `==`.

Finally, what I needed most in my project was structures. It's similar to R's list or Python's nested Dict, where it can store various kinds of information.

#### Strong Typing with Flexible Handling----Interface(interface)

Interfaces are a new concept I encountered when接触Go. It is quite different from the API I've heard before. I spent quite some time trying to understand it... here's a record.

First, according to the previous section, we know that Go is strongly typed. This strong typing runs through the entire language system and can sometimes cause unnecessary troubles. For example, when writing functions, you need to define the types of the receiving and returning variables, input and return variables must be given types. Once the type is determined, the function can only receive/this type of variable, it cannot change. This is very annoying. We write programs, there are many times when we have some functions that need to receive some objects with unclear situations, then according to the situation of the object for further processing. But under the rules of strong typing, receiving things must be a determined type. Then you have to write two duplicate functions to receive different types of objects and do similar processing.

Take an actual example from `vcfgo`'s `Info().Get()` function. This function's purpose is to take out corresponding information in the INFO column of a vcf file and return it. Each INFO item may have multiple values, if there is only one value, then return the corresponding string, if there are multiple, then put them into a list with strings as elements and return it. This looks very simple, but under the rules of strong typing, it cannot be directly implemented. Because when defining functions, you must define the type of the returned variable. A variable's type defined as a string cannot store a list of strings, vice versa. Then to achieve the target function, do we have to let the function have multiple return values? This is when `interface` comes into play. The definition of an interface is that as long as an object implements the methods defined in the interface, it has implemented this interface and can be passed in the form of this type. Returning back to the previous problem, just define an interface, then let the things you want to return implement this interface, so they can be returned in the form of interfaces.

Finally, note that since the return is an interface, their types are all `interface`. If you need to pass them further, you may still need to use assertion to take them out from the interface as the original type. For example, there's a line in my project:

```go
if annObj, ok := anns.(string); ok { // Single annotation, directly parse
	annTarget = annObj
}
```

This line of code where `ann` is the return value of the `Info().Get()` function. This interface may be a string or a list of strings, and I use `anns.(string)` to try to assert that it's a string inside. If the assertion is correct, then directly assign the parsed string, if incorrect, then perform another processing.

### File Operations
Python's file reading and writing are very simple. Just one line of code opens a file, and the obtained object is iterable. You can directly use `for` loop to process file content line by line. Writing files is similar, it can be done with just a few lines of code too. At the same time, Python has a `with` statement that automatically closes the file handle after operations are completed, which is very convenient.

```python
with open("file", "r") as f:
    for line in f:
        print(line)
```

Compared to Go, file reading and writing are much more troublesome. For example, if I want to read file content line by line and print it:

```go
file, err := os.Open(vcfFile)
defer file.Close() // Automatically close the file handle, similar to using 'with' in Python
if err != nil {
	panic(err)
}
reader := bufio.NewReader(file)
for {
	str, err := reader.ReadString('\n')
	if err == io.EOF {
		break
	}
	fmt.Println(str)
}
```

The main steps are first opening a file, then specifying with another function to read one line of content using `'\n'` as the delimiter. Then print the read content to the screen. Actually, the core part is similar to Python, but some things have become self-handled, so the code naturally becomes longer. For example, after opening the file, error handling in Go has to be done manually. If Python can't find a file when reading, it will automatically throw an exception and tell you that the file doesn't exist, then exit. In contrast, writing Go requires manual `panic` and displaying the error message from `err`. In Python, using `for` to iterate over the opened file object allows you to read file content line by line, and it will automatically terminate when reaching the end of the file. In Go, this process has to be written manually, even specifying that files are separated by `'\n'`.

Writing files is similar to reading, but slightly more complex in Go. Here's a pit I'll fill later.

### Nested Data Structures

Python commonly used data structures can be nested to form complex data structures to meet actual data processing needs. This kind of nested data structure is very flexible because dict/list/tuple/set have no special restrictions on what they store, so we can easily compose a series of different things into one organized structure. This is one reason why I prefer Python over R for certain data processing tasks (R has some internal structures that only store similar types of elements). For example, in the refGene file information stored in my project, my program needs to read out a gene's transcript number, exon start and end positions, transcription direction, gene name, etc. series of information. Among them, the transcript number is a string, the exon start and end positions store a list of integers, the transcription direction is a single character, and the gene name is also a string. This series of information can be stored in a large dictionary with `Key: Value` form, and you can get the corresponding value according to the Key when needed.

This is impossible in Go, because `list` corresponds to arrays, `dict` corresponds to maps, in Go they can only store things declared beforehand:

```go
outStrList := []string{"out1", "out2"} // Can only store strings
gene2tran := make(map[string][]string) // Key and Value must be strings
```

Therefore, to achieve similar nested data structures as Python, you have to use the structure `struct`:

```go
// Structure I defined in my project to store transcript information
type geneRecord struct {
	exon        [][]int64
	intron      [][]int64
	transRegion [2]int64
	cdsRegion   [2]int64
	chr         string
	strand      string
	geneSymbol  string
}
tran2info := make(map[string]geneRecord) // map of structures
```

### Control Structures

Go's control structures are even simpler than Python... judgment uses `if/else`, loop uses `for`, multi-branch structure uses `switch`------only these.

```go
// if else and for examples
if cols[3] == "+" {
	for idx, _ := range exonStarts {
		exons = append(exons, []int64{exonStarts[idx], exonEnds[idx], exonFrames[idx]})
		if idx != len(exonStarts)-1 {
			introns = append(introns, []int64{exonEnds[idx], exonStarts[idx+1]})
		}
	}
} else {
	for idx, _ := range exonStarts {
		exons = reverSlice(append(exons, []int64{exonStarts[idx], exonEnds[idx], exonFrames[idx]}))
		if idx != len(exonStarts)-1 {
			introns = reverSlice(append(introns, []int64{exonEnds[idx], exonStarts[idx+1]}))
		}
	}
}

// switch example
switch {
case AinA && BinB: // Target fusion, corresponding correct
	locates["breakA"] = locAinA
	locates["breakB"] = locBinB
	return locates, "right", nil
case AinB && BinA: // Target fusion, corresponding reversed
	locates["breakA"] = locAinB
	locates["breakB"] = locBinA
	return locates, "wrong", nil
default:
	err := errors.New("Not Fusion")
}
```

Although there are fewer things, they are enough...

Another thing to note is that Go provides `goto` as a jumping control statement besides `break` and `continue`. However, I haven't used it yet...

### Building Functions

Creating functions in Go is not special... just pay attention to declaring the types of input parameters and return values first (in the example below, `reverSlice` follows with the passed-in parameter, outside the parentheses are the declarations of the returned value). Once declared, they cannot be changed. If you need flexible inputs and outputs, you can rely on `interface`.

```go
func reverSlice(a [][]int64) [][]int64 {
	for i := len(a)/2 - 1; i >= 0; i-- {
		opp := len(a) - 1 - i
		a[i], a[opp] = a[opp], a[i]
	}
	return a
}
```

### Error Handling

Essentially, Go does not have Python's specialized error handling or exception throwing operation. However, when using Go, you will find that many functions return two values, for example, the file opening function:

```go
file, err := os.Open(yourFile)
```

The above file opening function returns not only the opened file object but also another variable. Personally, I understand this as something similar to the status code returned by executing a Shell command. Generally speaking, if the program execution within the function is normal, this variable will be `nil`. If there is an error during the function's part of execution, it will return an `error` class object containing the error description string. Processing this variable is equivalent to using `Try: Except:` in Python.

### String Handling

Go's string handling mainly relies on the `strings` package, which includes basic operations such as connection, splitting, replacement, and trimming of whitespace... the actual usage experience is similar to R's `stringr`... there are many things... but it feels not so convenient.

Another thing to note is that Python's formatted strings in Go do not have an equivalent. So you cannot create a large text template in Go and fill in the content.

### JSON Handling

Go also has a specialized `json` package to read JSON file contents into Go data structures, and convert Go objects into JSON strings. However, due to type restrictions, it may be slightly more complex than Python...

```go
tarTranFile := "RXA_gene_trans.json"
data, err := ioutil.ReadFile(tarTranFile)
if err != nil {
	return
}
var targetTransFile interface{}
err = json.Unmarshal(data, &targetTransFile)
targetTrans := targetTransFile.(map[string]interface{})
```

You can see that to read the written JSON file, I executed opening the file, loading the content, using a specific function to process it, and then because the return is `interface`, so finally did type assertion.

Similar operations... in Python, 1~2 lines can be completed...

### Program Compilation

Unlike Python which translates and executes code line by line, Go language needs to be compiled before use (although there's also `go run`), so you need to study how to compile executable files.
(Placeholder for digging a hole)

### Unique Features of Python and Go

This section records some unique features of Python and Go respectively (digging a hole placeholder).

#### Python Part

- List/Dictionary Comprehensions: In Python, list/dictionary comprehensions are very convenient iterable object generation/handling methods that can effectively simplify code (whether there's an efficiency boost isn't clear), but in Go, you have to write loops obediently.
- String Formatting: The benefit of Python's string formatting is that it allows you to fill in a large amount of text one by one. This writing method makes the final appearance of the string very clear and easy to understand. In contrast, Go can only concatenate strings.

#### Go Part

- Interface (see above)
- `switch` control structure
- Convenient parallelization: Although Python has multi-threading, multi-processing, coroutines... but they are not as efficient or convenient as Go. This is one of the reasons why I decided to learn Go.
```