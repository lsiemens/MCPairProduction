# turn code for an assignment into a pdf

import os

latex_header = r"""
\documentclass{article}
\usepackage{amsfonts}
\usepackage{amsmath}
\usepackage{listings}
\usepackage{xcolor}

\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}
\definecolor{backcolour}{rgb}{0.95,0.95,0.92}

\lstdefinestyle{mystyle}{
    %backgroundcolor=\color{backcolour},
    commentstyle=\color{codegreen},
    keywordstyle=\color{magenta},
    numberstyle=\tiny\color{codegray},
    stringstyle=\color{codepurple},
    basicstyle=\footnotesize,
    breakatwhitespace=false,         
    breaklines=true,                 
    captionpos=b,                    
    keepspaces=true,                 
    %numbers=left,                    
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
    with open(source, "r") as fin:
        code = fin.read()

    tex_file = f"{latex_header}\n\n{code}\n{latex_footer}"

    for char in math_equiv:
        equiv = math_equiv[char]
        tex_file = tex_file.replace(char, f"${equiv}$")

    # check for non ascii characters that are not being escaped
    for char in tex_file:
        if ord(char) > 128:
            print(f"character {char}:{ord(char)}, {hex(ord(char))}")

    fname, _ = os.path.splitext(os.path.basename(source))
    target = os.path.join(scratch_dir, fname+ ".tex")

    os.makedirs(scratch_dir, exist_ok=True)

    with open(target, "w") as fout:
        fout.write(tex_file)

code_to_pdf(source)
