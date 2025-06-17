Below is a two‑part “kit” you can copy‑and‑paste into your workflow:

Part	Purpose
A. Ruleset (for the LLM to follow)	Exhaustive checklist of CRAN‑compliant roxygen2 requirements, including the tricky S3 corner‑cases.
B. Prompt Skeleton (what you actually send)	Pre‑filled template you complete with your code‑snippet(s); it embeds the ruleset via a single instruction so the model always sees it.


⸻

A. Ruleset — CRAN‑compliant roxygen2 documentation

Strictly obey every numbered rule.
If any rule conflicts with the code or with another rule, follow the stricter rule.

1  General style

1.1 Use roxygen2 comment prefix #'  for every line.
1.2 Lines ≤ 80 characters; wrap long sentences.
1.3 Sections order: @title → @description → @param → @return → @details (opt.) → @examples → other tags.
1.4 Each exported object must have exactly one @export (or @exportS3Method …) tag.
1.5 Avoid \\dontrun{} altogether; use \\donttest{} if example is slow, external I/O, or uses internet.
1.6 If an object is not exported add @keywords internal or omit @export.

2  Content requirements

2.1  @title: ≤ 60 characters, sentence‑case, ends without period.
2.2  @description: 1–2 short sentences; no code examples.
2.3  @param: one entry per formal argument, in order, "arg" then concise explanation.
2.4  @return: → type (e.g. "data.frame"), shape/dimensions if relevant, plus 1‑sentence description.
2.5  @examples: runnable in clean R session (R CMD check --as‑cran), ≤ 5 sec CPU; wrap long-running lines in \\donttest{}.
2.6  Include @section Warning: only when a real hazard exists (CRAN discourages filler).
2.7  Use @seealso, @references, @family, etc. sparingly—only if they add navigational value.

3  S3 generics and methods

3.1  Generic defined in this package
  • Document generic and at least one method in one roxygen block.
  • Add @rdname <generic> to every method block to merge into a single .Rd.
  • @export on the generic, @method <generic> <class> + @export on each method.

3.2  Generic imported from another package (e.g. print)
  • Do not re‑document the generic.
  • Each method gets its own roxygen block or shares with sibling methods via @rdname, but must contain:

#' @method print myclass
#' @export

  • If the generic comes from a Suggested package, replace @export with
@exportS3Method pkg::print.

3.3  Use @inheritParams <generic> or @inherit to avoid repeating identical param docs across methods.

4  Namespace and imports (roxygen handles automatically if you add tags)

4.1  For S3 methods on generics in Imports, write @importFrom pkg generic in one block per file.
4.2  If @exportS3Method is used, no additional NAMESPACE directive required.

5  Data objects

5.1  Document every data/ object with a roxygen block in R/.
5.2  Use @format to describe columns / list components; end each bullet with a period.

6  Examples of tricky cases

(a) S3 generic you own)

#' Compute something amazing
#'
#' @description Generic for amazing computations.
#' @param x An object.
#' @param ... Passed to methods.
#' @return A numeric scalar.
#' @export
amaze <- function(x, ...) UseMethod("amaze")

(b) Method shares docs with generic)

#' @rdname amaze
#' @method amaze data.frame
#' @export
amaze.data.frame <- function(x, ...) { ... }

(c) Method for imported generic (e.g. print)

#' Print a \\code{myclass} object
#'
#' @description Pretty printer for objects of class \\code{myclass}.
#' @inheritParams base::print
#' @return Invisibly returns \\code{x}.
#' @method print myclass
#' @export
print.myclass <- function(x, ...) { ... }


⸻

B. Prompt Skeleton — copy, fill, send

You are roxygen‑bot 🦜.  
Your task: generate **roxygen2 documentation** that satisfies all rules A.1–A.6 below for the R code I will supply.

### Rules A (CRAN compliance)
<PASTE the numbered ruleset above verbatim or store in system‑prompt>

### My package context
• Package name:   <pkgName>
• Imported pkgs:  <imports, suggests>
• Objects to document (all exported):
  - <function / S3 generic / data etc.>
  - …

### Source code to document
```r
<PASTE R FUNCTIONS / DATA DEFINITIONS HERE>

Output format

Return only roxygen2 blocks (no narrative, no fencing), each placed immediately above the object it documents.
If multiple objects share documentation, use @rdname accordingly.
Do not emit .Rd code or NAMESPACE directives.

---

### How to use it

1. **Copy** the Ruleset once into your system‑prompt (or prepend it each time, if your LLM host lacks system‑prompt support).  
2. **Fill** the skeleton with package name, imports, and the raw R code.  
3. **Send** the whole prompt to the LLM.  
4. **Run** `devtools::document()` and `R CMD check --as‑cran`; zero NOTES/ERRORS should relate to documentation.

### Why this works

* The ruleset encodes every CRAN pain‑point (return tags, S3 edge‑cases, no `\\dontrun{}`, etc.).  
* The numbered format lets the model reference constraints unambiguously.  
* The output‑only instruction prevents chatty prose that would fail `R CMD check`.

---

**Citations**

Key guidance distilled from:

* roxygen2 official vignette on documentation ordering and S3 handling  [oai_citation:0‡cran.r-project.org](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html?utm_source=chatgpt.com)  
* roxygen2 namespace vignette for `@exportS3Method` and method registration  [oai_citation:1‡cran.r-project.org](https://cran.r-project.org/web/packages/roxygen2/vignettes/namespace.html?utm_source=chatgpt.com)  
* roxygen2 7.3.1 changelog confirming current syntax (2024‑01‑22)  [oai_citation:2‡roxygen2.r-lib.org](https://roxygen2.r-lib.org/news/index.html?utm_source=chatgpt.com)

Use this kit as‑is or adapt it to your house style—the checklist already meets CRAN’s 2025 policies.