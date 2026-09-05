%% wlpeerj.cls - Document class for PeerJ and Nature Portfolio npj submissions
%% Compatible with Overleaf / PeerJ template
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{wlpeerj}[2022/01/01 v1.3 PeerJ and npj submission class]

\DeclareOption{lineno}{\def\@lineno{true}}
\DeclareOption{fleqn}{\PassOptionsToPackage{fleqn}{amsmath}}
\DeclareOption{10pt}{\PassOptionsToClass{10pt}{article}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax

\LoadClass[10pt,twocolumn]{article}

\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{amsmath,amssymb}
\RequirePackage{graphicx}
\RequirePackage{booktabs}
\RequirePackage{xcolor}
\RequirePackage{geometry}
\RequirePackage{fancyhdr}
\RequirePackage{authblk}
\RequirePackage{microtype}
\RequirePackage{titlesec}
\RequirePackage{caption}
\RequirePackage[colorlinks=true, allcolors=blue]{hyperref}

% Geometry
\geometry{
  a4paper,
  top=2.0cm,
  bottom=2.0cm,
  left=1.5cm,
  right=1.5cm,
  columnsep=14pt,
  headheight=14pt
}

% Colors
\definecolor{headerblue}{RGB}{0, 86, 150}
\definecolor{headerrpm}{RGB}{199, 0, 57}
\definecolor{boxbg}{RGB}{244, 248, 251}
\definecolor{boxborder}{RGB}{208, 225, 253}

% Line numbering
\ifdefined\@lineno
  \RequirePackage{lineno}
  \AtBeginDocument{\linenumbers}
\fi

% Running Headers and Footers
\pagestyle{fancy}
\fancyhf{}
\rhead{\small\sffamily\thepage}
\lhead{\small\sffamily Manuscript in preparation --- target journal: npj Microgravity}
\rfoot{\footnotesize\sffamily\today}
\lfoot{\footnotesize\sffamily Barker et al. | Simulated Microgravity Adaptation in \textit{M. marinum}}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0.4pt}

% Section Titles Styling
\titleformat{\section}
  {\large\sffamily\bfseries\color{headerblue}}{\thesection}{1em}{}
\titleformat{\subsection}
  {\normalsize\sffamily\bfseries}{\thesubsection}{1em}{}
\titleformat{\subsubsection}
  {\small\sffamily\bfseries}{\thesubsubsection}{1em}{}

% Captions
\captionsetup{
  font={small,sf},
  labelfont={bf},
  labelsep=period
}

% Corresponding Author Macro
\newcommand{\corrauthor}[3][*]{%
  \def\@corrmark{#1}%
  \def\@corrname{#2}%
  \def\@corremail{#3}%
}

% Store abstract when defined in preamble.
% \gdef takes a parameter text then a braced body, so "\gdef\@abstract\bgroup"
% never stored anything - it made \bgroup the parameter text and let the abstract
% body be typeset immediately, which errors with "Missing \begin{document}".
% \RenewEnviron captures the environment body in \BODY without typesetting it.
\RequirePackage{environ}
\def\@abstract{}
\RenewEnviron{abstract}{\global\let\@abstract\BODY}

% Redefine maketitle to print title banner, authors, and abstract box in twocolumn
\renewcommand{\maketitle}{%
  \twocolumn[%
    \vspace{-0.4cm}%
    \noindent%
    \begin{minipage}{\textwidth}%
      {\LARGE\sffamily\bfseries\color{headerblue} \@title \par}%
      \vspace{0.3cm}%
      {\large\sffamily \@author \par}%
      \vspace{0.15cm}%
      \ifdefined\@corrname
        {\footnotesize\sffamily $^{*}$Correspondence: \@corrname\ifdefined\@corremail\ifx\@corremail\@empty\else\ (\href{mailto:\@corremail}{\@corremail})\fi\fi \par}%
      \fi
      \vspace{0.35cm}%
      \ifx\@abstract\@empty\else
        \noindent%
        \fcolorbox{boxborder}{boxbg}{%
          \begin{minipage}{\dimexpr\textwidth-2\fboxsep-2\fboxrule\relax}
            \vspace{0.12cm}%
            {\small\sffamily\bfseries\color{headerblue} ABSTRACT}\par\vspace{0.08cm}%
            {\small\itshape \@abstract \par}%
            \vspace{0.12cm}%
          \end{minipage}%
        }%
        \vspace{0.5cm}%
      \fi
    \end{minipage}%
  ]%
}

\endinput
