%% geoStat.cls
%% Based on ElegantBook, refined for GeoStatistics Reports
%% 2026/01/23 v4.0

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{geoStat}[2026/01/23 v4.0 GeoStatistics Book Class]

\RequirePackage{kvoptions}
\RequirePackage{etoolbox}

% ----- Options Setup -----
\SetupKeyvalOptions{family=GEO, prefix=GEO@, setkeys=\kvsetkeys}

\DeclareStringOption[blue]{color}
\DeclareStringOption[fancy]{mode}
\DeclareStringOption[newtx]{math}

% Color Presets
\DeclareVoidOption{green}{\kvsetkeys{GEO}{color=green}}
\DeclareVoidOption{cyan}{\kvsetkeys{GEO}{color=cyan}}
\DeclareVoidOption{blue}{\kvsetkeys{GEO}{color=blue}}
\DeclareVoidOption{gray}{\kvsetkeys{GEO}{color=gray}}
\DeclareVoidOption{black}{\kvsetkeys{GEO}{color=black}}

% Mode Presets
\DeclareVoidOption{fancy}{\kvsetkeys{GEO}{mode=fancy}}
\DeclareVoidOption{simple}{\kvsetkeys{GEO}{mode=simple}}

% Pass unknown options to ctexbook
\DeclareDefaultOption{\PassOptionsToClass{\CurrentOption}{ctexbook}}

\ProcessKeyvalOptions*\relax

% ----- Load Base Class -----
% zihao=-4: set font size to Xiaosi (approx 12pt)
% openany: chapters can start on any page
\LoadClass[a4paper,oneside,openany,UTF8,zihao=-4,fontset=none]{ctexbook}

% ----- Core Packages -----
\RequirePackage{amsmath,amssymb,amsfonts,mathrsfs}
\RequirePackage{bm}
\RequirePackage{graphicx}
\RequirePackage[export]{adjustbox}
\RequirePackage{float}
\RequirePackage{booktabs}
\RequirePackage{array}
\RequirePackage{makecell}
\RequirePackage{multicol,multirow}
\RequirePackage{siunitx}
\RequirePackage{xcolor}
\RequirePackage{fancyvrb}
\RequirePackage{hyperref}
\RequirePackage{geometry}
\RequirePackage{caption}
\RequirePackage{titlesec}
\RequirePackage{fancyhdr}
\RequirePackage{setspace}
\RequirePackage[many]{tcolorbox}
\RequirePackage{anyfontsize}
\RequirePackage{tikz}
\RequirePackage[shortlabels,inline]{enumitem}
\RequirePackage{subcaption}
\usetikzlibrary{backgrounds,calc,shadows,positioning,fit,arrows.meta,shapes.geometric}

\setlength{\emergencystretch}{1em} % 解决部分行溢出问题

% ----- Geometry -----
\geometry{
  a4paper,
  top=25.4mm, bottom=25.4mm,
  left=25mm, right=25mm,
  headheight=2.17cm,
  headsep=4mm,
  footskip=12mm
}

% ----- Fonts -----
\ifdefstring{\GEO@math}{newtx}{
  \RequirePackage{newtxmath}
  \RequirePackage{esint}
}{\relax}

% 强制显式设置中英文字体（适用于 XeLaTeX）
\RequirePackage[no-math]{fontspec}

% 英文采用 TeX Gyre Termes (Times New Roman 的高质量替代，通常随 TeX Live 安装)
\setmainfont{TeXGyreTermesX}[
    UprightFont = *-Regular,
    BoldFont = *-Bold,
    ItalicFont = *-Italic,
    BoldItalicFont = *-BoldItalic,
    Extension = .otf
]

\setCJKmainfont[
    Path = fonts/,
    Extension = .ttf,
    BoldFont = SimHei,
    ItalicFont = SimKai,
    BoldItalicFont = SimHei,
    BoldItalicFeatures = {FakeSlant=0.2},
    AutoFakeBold = true
]{SimSun}

\setCJKsansfont[Path = fonts/, Extension = .ttf]{SimHei}
\setCJKmonofont[Path = fonts/, Extension = .ttf]{SimFang}

% ----- Colors Definitions (From ElegantBook) -----
\definecolor{winered}{rgb}{0.5,0,0}

\ifdefstring{\GEO@color}{green}{
  \definecolor{structurecolor}{RGB}{0,120,2}
  \definecolor{main}{RGB}{0,120,2}
  \definecolor{second}{RGB}{230,90,7}
  \definecolor{third}{RGB}{0,160,152}
}{\relax}
\ifdefstring{\GEO@color}{cyan}{
  \definecolor{structurecolor}{RGB}{31,186,190}
  \definecolor{main}{RGB}{59,180,5}
  \definecolor{second}{RGB}{175,153,8}
  \definecolor{third}{RGB}{244,105,102}
}{\relax}
\ifdefstring{\GEO@color}{blue}{
  \definecolor{structurecolor}{RGB}{60,113,183}
  \definecolor{main}{RGB}{0,166,82}
  \definecolor{second}{RGB}{255,134,24}
  \definecolor{third}{RGB}{0,174,247}
}{\relax}
\ifdefstring{\GEO@color}{gray}{
  \definecolor{structurecolor}{RGB}{150,150,150}
  \definecolor{main}{RGB}{150,150,150}
  \definecolor{second}{RGB}{150,150,150}
  \definecolor{third}{RGB}{150,150,150}
}{\relax}
\ifdefstring{\GEO@color}{black}{
  \definecolor{structurecolor}{RGB}{0,0,0}
  \definecolor{main}{RGB}{0,0,0}
  \definecolor{second}{RGB}{0,0,0}
  \definecolor{third}{RGB}{0,0,0}
}{\relax}

% Cover line color usually matches 'second' color
\colorlet{coverlinecolor}{second}

% ----- Hyperref -----
\hypersetup{
  breaklinks,
  unicode,
  linktoc=all,
  bookmarksnumbered=true,
  bookmarksopen=true,
  colorlinks,
  linkcolor=winered,
  citecolor=winered,
  urlcolor=winered,
  plainpages=false,
  pdfstartview=FitH
}

% ----- 右上角水印 Logo 控制 -----
\newif\if@watermark\@watermarkfalse
\newcommand{\watermarklogo}[1]{\gdef\@watermarklogo{#1}\@watermarktrue}
\newcommand{\watermarkoff}{\@watermarkfalse}
\newcommand{\watermarkon}{\@watermarktrue}
\newcommand{\watermarksize}[1]{\gdef\@watermarksize{#1}}
\newcommand{\watermarkopacity}[1]{\gdef\@watermarkopacity{#1}}
\newcommand{\watermarkpos}[2]{\gdef\@watermarkx{#1}\gdef\@watermarky{#2}}
\gdef\@watermarksize{2.5cm}
\gdef\@watermarkopacity{0.8}
\gdef\@watermarkx{0cm}
\gdef\@watermarky{0cm}
\gdef\@watermarklogo{picture/logo/LongRuan1.jpg}

\newcommand{\@drawwatermark}{%
  \if@watermark%
    \begin{tikzpicture}[remember picture, overlay]
      % 基准点设在页眉右侧边界，且底部对齐页眉线
      % 页眉线位置：距离顶端为 top - headsep = 25.4mm - 4mm = 21.4mm
      \node[anchor=south east, xshift=\@watermarkx, yshift=\@watermarky, opacity=\@watermarkopacity] at ($(current page.north east) + (-25mm, -21.4mm)$) {
        \includegraphics[width=\@watermarksize, keepaspectratio]{\@watermarklogo}
      };
    \end{tikzpicture}%
  \fi%
}

% ----- Headers & Footers -----
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\color{structurecolor}\normalfont\leftmark}
\fancyhead[R]{\color{structurecolor}\normalfont\rightmark\@drawwatermark}
\fancyfoot[C]{\color{structurecolor}\small\thepage}
\renewcommand{\headrulewidth}{1pt}
\renewcommand{\headrule}{\color{structurecolor}\hrule width\textwidth}

% 覆盖 plain 样式以包含水印
\fancypagestyle{plain}{%
  \renewcommand{\headrulewidth}{0pt}%
  \fancyhf{}%
  \fancyhead[R]{\@drawwatermark}%
  \fancyfoot[C]{\color{structurecolor}\small\thepage}%
}

% 专门用于前言等无页眉文字但需水印的页面
\fancypagestyle{watermarkonly}{%
  \renewcommand{\headrulewidth}{0pt}%
  \fancyhf{}%
  \fancyhead[R]{\@drawwatermark}%
  \fancyfoot[C]{\color{structurecolor}\small\thepage}%
}

% ----- Sectioning Formatting -----
\ctexset{
  chapter = {
    format = \Huge\bfseries\color{structurecolor}\centering,
    name = {第,章},
    number = \arabic{chapter},
    beforeskip = -10pt,
    afterskip = 20pt,
  },
  section = {
    format = \Large\bfseries\color{structurecolor},
  },
  subsection = {
    format = \large\bfseries\color{structurecolor},
  }
}

% ----- Elegant Introduction Environment -----
\tcbset{
    introductionsty/.style={
        enhanced,
        breakable,
        colback=structurecolor!10,
        colframe=structurecolor,
        fonttitle=\bfseries,
        colbacktitle=structurecolor,
        fontupper=\itshape,
        attach boxed title to top center={yshift=-3mm,yshifttext=-1mm},
        boxrule=0pt,
        toprule=0.5pt,
        bottomrule=0.5pt,
        top=8pt,
        before skip=8pt,
        sharp corners
    },
}

% The square icon for the list
\newcommand*{\eitemi}{\tikz \draw [baseline, ball color=structurecolor,draw=none] circle (2pt);}

\newenvironment{introduction}[1][内容提要]{
  \begin{tcolorbox}[introductionsty,title={#1}]
    \begin{multicols}{2}
      \begin{itemize}[label={\eitemi}, nosep]
}{
      \end{itemize}
    \end{multicols}
  \end{tcolorbox}
}

% ----- Theorem Environments (Tcolorbox) -----
\tcbset{
    common/.style={
        lower separated=false,
        coltitle=white,
        colback=gray!5,
        boxrule=0.5pt,
        fonttitle=\bfseries,
        enhanced,
        breakable,
        top=8pt,
        before skip=8pt,
        attach boxed title to top left={yshift=-0.11in, xshift=0.15in},
        boxed title style={boxrule=0pt, colframe=white, arc=0pt, outer arc=0pt}
    },
    defstyle/.style={common, colframe=main, colback=main!5, colbacktitle=main},
    thmstyle/.style={common, colframe=second, colback=second!5, colbacktitle=second},
    prostyle/.style={common, colframe=third, colback=third!5, colbacktitle=third},
    notestyle/.style={common, colframe=gray!60!black, colback=gray!5, colbacktitle=gray!60!black}
}

% Define theorems
\newtcbtheorem[number within=chapter]{theorem}{定理}{thmstyle}{thm}
\newtcbtheorem[use counter from=theorem]{lemma}{引理}{thmstyle}{lem}
\newtcbtheorem[use counter from=theorem]{definition}{定义}{defstyle}{def}
\newtcbtheorem[use counter from=theorem]{example}{例}{prostyle}{exa}
\newtcbtheorem[use counter from=theorem]{note}{注}{notestyle}{not}
\newtcbtheorem[use counter from=theorem]{warning}{警告}{thmstyle}{warn}

% ----- Bibliography -----
\RequirePackage[backend=biber, style=gb7714-2015, maxcitenames=1, mincitenames=1]{biblatex}

% 配置 biblatex 支持带中文前缀的引用（针对 \citep 等命令）
% 使 postnote 直接显示原始输入（不自动添加 'p.' 或 'pp.'）
\DeclareFieldFormat{postnote}{#1}
\DeclareFieldFormat{multipostnote}{#1}

% 重新定义 \citep 命令以支持中文前缀
\DeclareCiteCommand{\citep}[\mkbibparens]
  {\usebibmacro{prenote}}
  {\usebibmacro{citeindex}%
   \bibhyperref{%
     \ifnameundef{labelname}
       {\printfield{labeltitle}}
       {\printnames{labelname}%
        \ifnumgreater{\value{labelname}}{1}
          {\space\iffieldequalstr{langid}{chinese}{等}{et~al.}}
          {}}%
     \setunit{\addcomma\space}%
     \printfield{labelyear}}}
  {\multicitedelim}
  {\iffieldundef{postnote}
    {}
    {\setunit{\addcomma\space}%
     \printfield{postnote}}}

% ----- Algorithms -----
\RequirePackage[ruled,linesnumbered,vlined,algochapter]{algorithm2e}
\SetAlgorithmName{算法}{算法}{算法列表}
\SetAlgoCaptionSeparator{:\ }
\SetAlCapHSkip{0pt}

% ----- Caption Settings -----
\captionsetup{
  font={small},
  labelfont={bf,color=structurecolor},
  labelsep=quad,
  figurename={图},
  tablename={表}
}

% ----- Cover Page -----
\newcommand{\subtitle}[1]{\gdef\@subtitle{#1}}
\newcommand{\institute}[1]{\gdef\@institute{#1}}
\newcommand{\version}[1]{\gdef\@version{#1}}
\newcommand{\cover}[1]{\gdef\@cover{#1}}
\newcommand{\logo}[1]{\gdef\@logo{#1}}
\newcommand{\logosize}[1]{\gdef\@logosize{#1}}
\newcommand{\secondlogo}[1]{\gdef\@secondlogo{#1}}
\newcommand{\secondlogosize}[1]{\gdef\@secondlogosize{#1}}

\renewcommand*{\maketitle}{%
  \hypersetup{pageanchor=false}
  \begin{titlepage}
    \newgeometry{margin=0in}
    \parindent=0pt
    % 1. 顶部封面图或色块
    {\setlength{\fboxsep}{0pt}%
    \ifcsname @cover\endcsname
        \begin{tikzpicture}
            \clip (0,0) rectangle (\linewidth, 0.45\textheight);
            \node[inner sep=0pt, anchor=center] at (0.5\linewidth, 0.225\textheight) {%
                \includegraphics[min width=\linewidth, min height=0.45\textheight, keepaspectratio]{\@cover}%
            };
        \end{tikzpicture}
    \else
        \colorbox{structurecolor!15}{\makebox[\linewidth][c]{\rule{0pt}{0.45\textheight}}}
    \fi
    
    % 2. 结构装饰条
    \colorbox{coverlinecolor}{\makebox[\linewidth][c]{\rule{0pt}{2mm}}}%
    }
    
    \vfill
    
    % 3. 标题与副标题区
    \hspace{2.5em}
    \begin{minipage}{0.85\linewidth}
        {\Huge\bfseries\color{structurecolor} \@title \par}
        \vspace{0.6cm}
        {\Large\bfseries\color{gray!80!black} \ifcsname @subtitle\endcsname \@subtitle \fi \par}
    \end{minipage}
    
    \vspace{2.5cm}
    
    % 4. 个人/项目信息区
    \hspace{3.5em}
    \begin{minipage}[t]{0.6\linewidth}
        \renewcommand{\arraystretch}{1.5}
        \color{gray!60!black}
        \begin{tabular}{l}
            \ifx\@author\empty\else {\Large \textbf{作者：} \@author} \\ \fi
            \ifcsname @institute\endcsname {\Large \textbf{单位：} \@institute} \\ \fi
            \ifx\@date\empty\else {\Large \textbf{日期：} \@date} \\ \fi
            \ifcsname @version\endcsname {\Large \textbf{版本：} \@version} \\ \fi
        \end{tabular}
    \end{minipage}
    
    % 5. Logo 放置 (右下角绝对定位)
    \begin{tikzpicture}[remember picture,overlay]
        \node[anchor=south east, xshift=-1.2in, yshift=1.2in, align=center] at (current page.south east) {
            \ifcsname @logo\endcsname
                \includegraphics[width=\ifcsname @logosize\endcsname\@logosize\else 4.5cm\fi, keepaspectratio]{\@logo}
            \fi
            \ifcsname @secondlogo\endcsname
                \\[0.1cm]
                \includegraphics[width=\ifcsname @secondlogosize\endcsname\@secondlogosize\else 4.5cm\fi, 
                keepaspectratio]{\@secondlogo}
            \fi
        };
    \end{tikzpicture}
    
    \vfill
    
    % 6. 底部装饰线
    \begin{center}
       \color{structurecolor}\rule{0.9\linewidth}{1.5pt}
    \end{center}
    \vspace{1.2cm}
    
    \restoregeometry
    \thispagestyle{empty}
  \end{titlepage}
  \hypersetup{pageanchor=true}
  \newpage
}

\endinput
%% End of geoStat.cls