# turn code for an assignment into a pdf

import os

latex_header = r"""
\documentclass{article}
\usepackage{amsfonts}
\usepackage{amsmath}
\usepackage{listings}
\usepackage{xcolor}

\definecolor{ccomment}{rgb}{0.5, 0.0, 0.2}
\definecolor{ckeyword}{rgb}{0.0, 0.2, 0.5}
\definecolor{cstring}{rgb}{0.0, 0.5, 0.2}

\lstdefinestyle{mystyle}{
    commentstyle=\color{ccomment},
    keywordstyle=\color{ckeyword},
    stringstyle=\color{cstring},
    breakatwhitespace=false,
    breaklines=true,
    captionpos=b,
    keepspaces=true,
    numbersep=5pt,
    showspaces=false,
    showstringspaces=false,
    showtabs=false,
    tabsize=2
}

\lstset{style=mystyle}
\begin{document}
\lstset{language=Python}

\begin{lstlisting}[mathescape=true]
"""

math_equiv = {"²":r"^2", "β":r"\beta", "σ":r"\sigma", "𝜋":r"\pi"}

latex_footer = "\end{lstlisting}\end{document}"

scratch_dir = "./scratch"

source = "./Assignment2.py"

def code_to_pdf(source):
    with open(source, "r", encoding="utf-8") as fin:
        code = fin.read()

    tex_file = f"{latex_header}\n\n{code}\n{latex_footer}"

    for char in math_equiv:
        equiv = math_equiv[char]
        tex_file = tex_file.replace(char, f"$\color{{ccomment}}{equiv}$")

    # check for non ascii characters that are not being escaped
    for char in tex_file:
        if ord(char) > 128:
            print(f"character {char}:{ord(char)}, {hex(ord(char))}")

    fname, _ = os.path.splitext(os.path.basename(source))
    target = os.path.join(scratch_dir, fname+ ".tex")

    os.makedirs(scratch_dir, exist_ok=True)

    with open(target, "w", encoding="utf-8") as fout:
        fout.write(tex_file)

code_to_pdf(source)
