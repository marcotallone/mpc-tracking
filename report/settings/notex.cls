%┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓%
%┃                                ~  NOTEX ~                                  ┃%
%┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{notex}[2024/05/17 A modern note-taking LaTeX class.]

\LoadClass[a4paper, 10pt]{article}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                 Geometry                                   │%
%└────────────────────────────────────────────────────────────────────────────┘%

\RequirePackage{geometry}

% Margins / Borders
\geometry{
  left    =   3.8cm,
  right   =   3.8cm,
  top     =   0.5cm,
  bottom  =   0.5cm,
  includeheadfoot
}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                  General                                   │%
%└────────────────────────────────────────────────────────────────────────────┘%

% Multiple columns
\RequirePackage{multicol}

% Lists (enumeration, ...)
\RequirePackage{enumitem}

% Hyperlinks
\RequirePackage{hyperref}

% Index
% \RequirePackage{index}

% Abstract
\RequirePackage{abstract}

% PDF inclusion
% \RequirePackage{pdfpages}

% Emoji support for LaTeX
% \RequirePackage{emoji}
% \setemojifont{Noto Color Emoji}

% Fontawesome5 package for icons
% \RequirePackage{fontawesome5}

% tcolorbox for colored boxes
% [breakable] option allows the box to be split across pages
% \RequirePackage[breakable]{tcolorbox}

% To suppress annoying overfull hbox warnings
\hbadness=10000

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                   Style                                    │%
%└────────────────────────────────────────────────────────────────────────────┘%

% No indentation
\setlength{\parindent}{0pt}

% NOTE: This is a quick set-up I use to print my name, the section and date
%       on all pages.
%       To use it uncomment it and define the commands:
%       \newcommand{\pastetitle}{MyTitle}
%       \newcommand{\pasteauthor}{MyName}
%       \newcommand{\pastedate}{\today}
%       in the preamble of your document, possibly before the \maketitle.

% Headers and footers
% \RequirePackage{fancyhdr}
% \pagestyle{fancy}
% \fancyhf{}
% \fancyhead[L]{\nouppercase{\footnotesize{\gray{\leftmark}}}}
% \fancyhead[C]{\nouppercase{\footnotesize{\gray{\pastetitle}}}}
% \fancyhead[R]{\nouppercase{\footnotesize{\gray{\pasteauthor\ | \pastedate}}}}
% \renewcommand{\headrulewidth}{0pt} % No horizontal line
% \fancyfoot[C]{\thepage} % Page number in the center

% Section titles
\RequirePackage{titlesec}
\titleformat*{\section}{\fontsize{18}{16}\selectfont\bfseries}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                   Fonts                                    │%
%└────────────────────────────────────────────────────────────────────────────┘%

% Use Helvetica font
% \RequirePackage[scaled]{helvet}

% Use the Tex Gyre Heros font
% \RequirePackage{tgheros}

% Use the Computer Modern Sans Serif font
% \RequirePackage{cmbright}

% Use the Latin Modern Sans Serif font (default)
% \RequirePackage{lmodern}

% Set the default font to sans-serif
% \renewcommand*\familydefault{\sfdefault}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                   Math                                     │%
%└────────────────────────────────────────────────────────────────────────────┘%

\RequirePackage{mathtools}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{amsfonts}

% To display nice fractions in text mode
\RequirePackage{nicefrac}

% To use \mathbbm{1} to display the identity matrix
\RequirePackage{bbm}

% To use \mathcal{O} to display big O notation
% \RequirePackage{mathrsfs}

% To represent simplifications with the \cancel{} command
\RequirePackage{cancel}

% To reset environments counters at every section
% (amsmath package required)
\numberwithin{figure}{section}
\numberwithin{equation}{section}
\numberwithin{table}{section}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                  Colors                                    │%
%└────────────────────────────────────────────────────────────────────────────┘%

\RequirePackage[table]{xcolor}
\RequirePackage{pagecolor}

% Theme

% Light Theme
\DeclareOption{light}{

  % Background and boxes
  \definecolor{background}{HTML}{FFFFFF}
  \definecolor{boxcolor}{HTML}{EFF1F3}
  \definecolor{darkred}{HTML}{D72323}
  \definecolor{darkblue}{HTML}{0A84FF}
  \definecolor{darkgreen}{HTML}{00C853}
  \definecolor{darkyellow}{HTML}{FFCB30}

  % Text color / shades
  \definecolor{black}{HTML}{000000}
  \definecolor{darkgray}{HTML}{161B22}
  \definecolor{gray}{HTML}{89929B}
  \definecolor{lightgray}{HTML}{C6CDD5}
  \definecolor{white}{HTML}{FFFFFF}
  \definecolor{textcolor}{HTML}{000000}       %(black)

  % Colors                                    % --HEX--
  \definecolor{red}{HTML}{BD2C00}             % #BD2C00
  \definecolor{orange}{HTML}{C9510C}          % #C9510C
  \definecolor{yellow}{HTML}{FFCB30}          % #FFCB30
  \definecolor{green}{HTML}{6CC644}           % #6CC644
  \definecolor{azure}{HTML}{54AEFF}           % #54AEFF
  \definecolor{blue}{HTML}{4078C0}            % #4078C0
  \definecolor{purple}{HTML}{6E5494}          % #6E5494

  % Background and text color
  \pagecolor{background}                      % Page Background
  \color{textcolor}                           % Main text color
  \definecolor{captionscolor}{HTML}{89929B}   % Caption colors
  \definecolor{iconscolor}{HTML}{89929B}      % Icons color
  \definecolor{linescolor}{HTML}{89929B}      % Lines color

}

% Dark Theme
\DeclareOption{dark}{

  % Background and boxes
  \definecolor{background}{HTML}{0D1117}
  \definecolor{boxcolor}{HTML}{21262D}
  \definecolor{darkred}{HTML}{D72323}
  \definecolor{darkblue}{HTML}{0A84FF}
  \definecolor{darkgreen}{HTML}{00C853}
  \definecolor{darkyellow}{HTML}{FFCB30}

  % Text shades
  \definecolor{black}{HTML}{000000}
  \definecolor{darkgray}{HTML}{161B22}
  \definecolor{gray}{HTML}{89929B}
  \definecolor{lightgray}{HTML}{C6CDD5}
  \definecolor{white}{HTML}{FFFFFF}
  \definecolor{textcolor}{HTML}{FFFFFF}       % (white)      

  % Colors                                    % --HEX--
  \definecolor{red}{HTML}{FA7970}             % #FA7970 
  \definecolor{orange}{HTML}{FAA356}          % #FAA356
  \definecolor{yellow}{HTML}{FFCB30}          % #FFCB30
  \definecolor{green}{HTML}{7CE38B}           % #7CE38B
  \definecolor{azure}{HTML}{A2D2FB}           % #A2D2FB
  \definecolor{blue}{HTML}{77BDFB}            % #77BDFB
  \definecolor{purple}{HTML}{CEA5FB}          % #CEA5FB

  % Background and text color
  \pagecolor{background}                      % Page Background
  \color{textcolor}                           % Main text color
  \definecolor{captionscolor}{HTML}{C6CDD5}   % Caption Colors
  \definecolor{iconscolor}{HTML}{FFFFFF}      % Icons color
  \definecolor{linescolor}{HTML}{C6CDD5}      % Lines color

}

% Set the default option
\ExecuteOptions{light}
\ProcessOptions\relax

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                  Figures                                   │%
%└────────────────────────────────────────────────────────────────────────────┘%

% General figures packages
\RequirePackage{graphicx}
\RequirePackage{wrapfig}

% Either one or the other commented!
% \RequirePackage{subfig}
\RequirePackage{subcaption}

\RequirePackage{adjustbox} % Scale figures
\RequirePackage[section]{placeins} % Figures inside sections

% Captions
\RequirePackage[labelfont={color=textcolor, bf},
            textfont={color=textcolor},
            margin=1cm]{caption}
\captionsetup{format=hang}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                  Tables                                    │%
%└────────────────────────────────────────────────────────────────────────────┘%

\RequirePackage{booktabs}
\RequirePackage{array}
\RequirePackage{multirow}
\RequirePackage{tabularx}

% Table templates (TT)
% ------------------------------------------------------------------------------
% Normal Table with tabular

% \begin{table}[htb]
%   \renewcommand{\arraystretch}{1.5} % Row height
%   \centering
%   \begin{tabular}{|c|c|c|}

%     % Header (different color)
%     \hline
%     \rowcolor{boxcolor}
%     \textbf{Column1} &
%     \textbf{Column2} &
%     \textbf{Column3} \\ 

%     % Rows
%     \hline
%     Row1 & 
%     Row1 & 
%     Row1 \\
%     \hline

%   \end{tabular}
%   \caption{MyTableCaption}
%   \label{tab:mytablelabel}
%   \renewcommand{\arraystretch}{1} % Reset row height to default
% \end{table}

% ------------------------------------------------------------------------------
% Table with tabularx

% \begin{table}[htb]
%   \renewcommand{\arraystretch}{1.5} % Row height
%   \centering
%   \begin{tabularx}{\textwidth}{|X|X|X|}

%     % Header (different color)
%     \hline
%     \rowcolor{boxcolor}
%     \textbf{Column1} &
%     \textbf{Column2} &
%     \textbf{Column3} \\ 

%     % Rows
%     \hline
%     Row1 & 
%     Row1 & 
%     Row1 \\
%     \hline

%   \end{tabularx}
%   \caption{MyTablexCaption}
%   \label{tab:mytablexlabel}
%   \renewcommand{\arraystretch}{1} % Reset row height to default
% \end{table}

% ------------------------------------------------------------------------------
% // TODO: maybe build environment/macros to deal with these tables...
% ------------------------------------------------------------------------------

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                                   TikZ                                     │%
%└────────────────────────────────────────────────────────────────────────────┘%

% TikZ packages
\RequirePackage{tikz}
\RequirePackage{scalerel}
\RequirePackage{pict2e}
\RequirePackage{tkz-euclide}
\RequirePackage{circuitikz} % To draw circuits

% Libraries
\usetikzlibrary{calc}
\usetikzlibrary{patterns,arrows.meta}
\usetikzlibrary{shadows}
\usetikzlibrary{external}

% PGFPlots
\RequirePackage{pgfplots}
\pgfplotsset{compat=newest}
\usepgfplotslibrary{statistics}
\usepgfplotslibrary{fillbetween}
% Standard, for all pgfplot:
\pgfplotsset{
	standard/.style={
		axis line style = thick,
		trig format=rad,
		enlargelimits,
		axis x line=middle,
		axis y line=middle,
		enlarge x limits=0.15,
		%enlarge y limits=0.15,
		every axis x label/.style={
      at={(current axis.right of origin)},
      anchor=north west
    },
		every axis y label/.style={
      at={(current axis.above origin)},
      anchor=south east
    },
		%grid=both
	}
}

% Memory from my thesis ❤
\definecolor{blueish}{RGB}{125,100,255}
\definecolor{yucky}{HTML}{808000}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                             Table of Contents                              │%
%└────────────────────────────────────────────────────────────────────────────┘%

% Set the depth of the table of contents to show:
% 0 - chapter
% 1 - section
% 2 - subsection
% 3 - subsubsection
\setcounter{tocdepth}{2}

% Set the depth of the sections numbering:
\setcounter{secnumdepth}{3}

% Table of contents
%\RequirePackage{tocloft}
\RequirePackage[subfigure]{tocloft}

% For filling dots between toc elements and page numbers
% \renewcommand{\cftpartleader}{\cftdotfill{\cftdotsep}}  % for parts
%\renewcommand{\cftchapleader}{\cftdotfill{\cftdotsep}}   % for chapters
% \renewcommand{\cftsecleader}{\cftdotfill{\cftdotsep}}   % for sections, if
                                                          % you really want!

% Remove page numbers for parts in the table of contents
\cftpagenumbersoff{part}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                             External Packages                              │%
%└────────────────────────────────────────────────────────────────────────────┘%

% Code
\RequirePackage{settings/code}

% Callouts 
% \RequirePackage{settings/callouts}

%┌────────────────────────────────────────────────────────────────────────────┐%
%│                             Macros / Commands                              │%
%└────────────────────────────────────────────────────────────────────────────┘%
% Personal commands to make life easier

% Quick \textbf{} bold font
\newcommand{\bft}[1]{\textbf{#1}}

% Quick \textit{} italic font
\newcommand{\itt}[1]{\textit{#1}}

% Quick \texttt{} monospace font
\newcommand{\ttt}[1]{\texttt{#1}}

% Quick \mathbf{} bold math font
\newcommand{\mbf}[1]{\mathbf{#1}}

% Quick \mathbb{} blackboard bold math font
\newcommand{\mbb}[1]{\mathbb{#1}}

% Quick \mathcal{} calligraphic math font
\newcommand{\mcal}[1]{\mathcal{#1}}

% Identity matrix
\newcommand{\id}{\mathbbm{1}}

% Quick colors
\newcommand{\red}[1]{\textcolor{red}{#1}}
\newcommand{\orange}[1]{\textcolor{orange}{#1}}
\newcommand{\yellow}[1]{\textcolor{yellow}{#1}}
\newcommand{\green}[1]{\textcolor{green}{#1}}
\newcommand{\azure}[1]{\textcolor{azure}{#1}}
\newcommand{\blue}[1]{\textcolor{blue}{#1}}
\newcommand{\purple}[1]{\textcolor{purple}{#1}}
\newcommand{\gray}[1]{\textcolor{gray}{#1}}

% Quick bold font colors
\newcommand{\bred}[1]{\textcolor{red}{\bft{#1}}}
\newcommand{\borange}[1]{\textcolor{orange}{\bft{#1}}}
\newcommand{\byellow}[1]{\textcolor{yellow}{\bft{#1}}}
\newcommand{\bgreen}[1]{\textcolor{green}{\bft{#1}}}
\newcommand{\bazure}[1]{\textcolor{azure}{\bft{#1}}}
\newcommand{\bblue}[1]{\textcolor{blue}{\bft{#1}}}
\newcommand{\bpurple}[1]{\textcolor{purple}{\bft{#1}}}
\newcommand{\bgray}[1]{\textcolor{gray}{\bft{#1}}}

% Command to print an horizontal line given the width in input (optional,
% default is \textwidth)
\newcommand{\oline}[1][\textwidth]{\noindent\rule{#1}{0.1pt}}
