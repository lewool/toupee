# Style Guidelines

The purpose of this document is to establish style guidelines to follow when adding and reviewing new files or MATLAB code to this repository.

[Richard Johnson](https://uk.mathworks.com/matlabcentral/profile/authors/22731-richard-johnson) writes, "Style guidelines are not commandments. Their goal is simply to help programmers write well." Well-written code implies code that is easy to read. Code that is easy to read is typically written in a consistent style, so new code should be as consistent as possible with the rest of the repository.

A MATLAB file's header documentation (aka a docstring) is written as follows:
* Functions have the following sections in the given order:
	- One-line summary description, typically as a verb phrase
	- Long description (optional)
	- Inputs (optional)
	- Outputs
	- Examples
	- See also (optional)
	- Warnings/Exceptions (optional)
	- Additional notes (optional)
	- Todos (optional)
* Classes have the following sections in the given order:
	- One-line summary description, typically as a noun phrase
	- Long description (optional)
	- Examples
	- See also (optional)
	- Warnings/Exceptions (optional)
	- Additional notes (optional)
	- Todos (optional)

A file's body documentation should adhere to the following:
* For classes, each method, property, and event should be documented. Method docstrings only need contain the 'One-line summary description' and 'Examples' sections. Additional explanations should be given for methods, properties, and events that are within blocks that are given non-default attributes. e.g.
```
properties (Dependent)
  % `Dependent` because...
  property1
end

methods (Access=protected)
  % `Access=protected` because...
  method1
end

events (ListenAccess=private)
  % `ListenAccess=private` because...
  event1
end
```
* Whitespace conventions:
	- A tab/indent is set at two spaces.
	- Whitespace between lines is used sparingly, but can be used to improve readability between blocks of code.
	- Each line contains no more than 75 characters. Whitespace should be used to align statements that span multiple lines via a hanging indent or single indent e.g.
	```
	stash = ['stash push -m "stash working changes before '...
             'scheduled git update"'];
	```
* Block quotes should be written as full sentences above the corresponding code. e.g.
  ```
  % Check that the inputs are properly defined.
  if len(varargin) ~= 1 && len(varargin) ~= 2 ...
  ```
* Inline quotes should be written as a short phrase that starts two spaces after the corresponding code. e.g.
  ```
  if len(varargin) == 1  % return input
  ```
* Variables should be documented where they are declared. e.g.
  ```
  % The Rigbox root directory
  root = getOr(dat.paths, 'rigbox');
  ```
* Variable names referenced in comments are surrounded by back ticks for readability e.g.
  ```
  % New signal carrying `a` with its rows flipped in the up-down direction
  b = flipud(a);
  ```
* In general, clarity > brevity. Don't be afraid to spread things out over a number of lines and to add block and in-line comments. Long variable names are often much clearer (e.g. `inputSensorPosCount` instead of `inpPosN`).
* ["Scare quotes"](https://www.chicagomanualofstyle.org/qanda/data/faq/topics/Punctuation/faq0014.html) are generally to be avoided.

Additional conventions:
* Naming conventions:
	- Variable and function names are in lower camelCase (e.g. `expRef`, `inferParameters`)
	- Class and class property names are in upper CamelCase (e.g. `AlyxPanel`, `obj.Token`)
* Same variable names used across files should have the same meanings, and vice versa.
* When assigning a function's output(s), use the name(s) defined in that function's header.

See http://www.datatool.com/downloads/MatlabStyle2%20book.pdf for further details on MATLAB style guidelines.