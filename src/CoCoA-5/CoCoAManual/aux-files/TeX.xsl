<?xml version="1.0" encoding="ascii"?>

<xsl:stylesheet version="1.1"
		xmlns:xsl="http://www.w3.org/1999/XSL/Transform">


<!--xslt document for creating the man directory for CoCoA-->

  <xsl:output method="text"/>
  <!-- <xsl:strip-space elements="*"/>	 -->
  <xsl:preserve-space elements="quotes"/>

<!-- ================================================================== -->

  <xsl:param name="CoCoAVersion" required="yes"/>

  <xsl:template match="help">\documentclass[a4paper]{mybook}
\ifx\pdfoutput\undefined % if we're not running pdftex
\else
\usepackage{hyperref}
\hypersetup{
  pdfpagelayout=OneColumn,
  pdfstartview      = FitH,%     % La page prend toute la largeur}
%  bookmarks, 
  bookmarksopenlevel = 1,
%  bookmarksnumbered,
%  pagecolor = {green},
  colorlinks = {true},
  linkcolor = {blue},
  naturalnames,
}
\fi
\usepackage[usenames,dvipsnames]{color}
\usepackage{fancyvrb} % coloured boxes around verbatim region
\usepackage{supertabular} % page breaks in supertabular
\usepackage{amsfonts}
\usepackage{longtable}

\def\cocoa{{CoCoA}}
\newenvironment{command}{}{} % does nothing: may be useful for debugging
\newcommand\SeeAlso{\par\textcolor{OrangeRed}{\textbf{\large See Also: }}}

\title{\Huge{\textsf{{\cocoa-}<xsl:value-of select="$CoCoAVersion"/>\\Manual}}}

\author{John Abbott and Anna M. Bigatti}


\setcounter{secnumdepth}{1} 
\setcounter{tocdepth}{1} 

\setlength{\textwidth}{17cm}
\setlength{\oddsidemargin}{-.5cm}
\setlength{\evensidemargin}{-.5cm}

\setlength{\textheight}{24cm}
\setlength{\topmargin}{-1.5cm}

\setlength{\parskip}{4pt}
\setlength{\parindent}{0pt}

\tolerance=5000

\renewcommand{\thechapter}{\thepart-\arabic{chapter}}

\begin{document}

\maketitle

\tableofcontents

<!-- ===============  MANUAL  =============== -->

<!-- ===============  COMMANDS  =============== -->

% -- NEW PART  --------------------------------
\part{Alphabetical List of Commands}
\setcounter{chapter}{-1}

<xsl:for-each select="cocoa_commands/chapter_letter">
\chapter{<xsl:apply-templates select="title"/>}  %----=== CHAPTER LETTER ===----
\label{<xsl:apply-templates select="title" mode="strip"/>}
<xsl:for-each select="command">
  <xsl:sort select="title"/>
\section{<xsl:apply-templates select="title"/>}
\label{<xsl:apply-templates select="title" mode="strip"/>}
\begin{command} % -- COMMAND: <xsl:apply-templates select="title"/> ------------
<!--\subsection*{Syntax}-->
<xsl:if test="syntax">
\begin{Verbatim}[label=syntax, rulecolor=\color{MidnightBlue},
frame=single]<xsl:apply-templates select="syntax"/>\end{Verbatim}
</xsl:if>

<!--\subsection*{Description}-->
<xsl:apply-templates select="description"/>
<xsl:if test="seealso/see">
\SeeAlso %---- SEE ALSO ----
  <xsl:for-each select="seealso/see">
    <xsl:apply-templates/>(\ref{<xsl:call-template name="title-strip"
    select="."/>} pg.\pageref{<xsl:call-template name="title-strip"
    select="."/>})<xsl:if test="not(position()=last())">, 
    </xsl:if>

  </xsl:for-each>

</xsl:if>
\end{command} % -- end command --------------------------------
</xsl:for-each>
</xsl:for-each>


<!-- enable/disable man about chapters -->
  <xsl:if test="true()">
% -- PART 1 --------------------------------

  <xsl:for-each select="manual_parts/part">
\part{<xsl:apply-templates select="title"/>}
\setcounter{chapter}{0}
    <xsl:variable name="part_num" select="position()"/>
    <xsl:for-each select="chapter">

% -- CHAPTER --------------------------------
\chapter{<xsl:apply-templates select="title"/>}
\label{<xsl:apply-templates select="title" mode="strip"/>}

      <xsl:variable name="chap_num" select="position()"/>
      <xsl:for-each select="section">

% -- SECTION --------------------------------
\section{<xsl:apply-templates select="title"/>}
\label{<xsl:apply-templates select="title" mode="strip"/>}

        <xsl:variable name="sect_num" select="position()"/>
        <xsl:apply-templates select="description"/>
<xsl:if test="seealso/see">
\SeeAlso %---- SEE ALSO ----
  <xsl:for-each select="seealso/see">
    <xsl:apply-templates/>(\ref{<xsl:call-template name="title-strip"
    select="."/>} pg.\pageref{<xsl:call-template name="title-strip"
    select="."/>})<xsl:if test="not(position()=last())">, 
    </xsl:if>

  </xsl:for-each>

</xsl:if>

      </xsl:for-each>
    </xsl:for-each>
  </xsl:for-each>

  </xsl:if>
<!-- enable/disable man about chapters -->

\end{document}
  </xsl:template>

<!-- ================================================================== -->

  <!-- taken from http://www.xml.com/pub/a/2002/06/05/transforming.html -->
  <xsl:template name="str_replace">
    <xsl:param name="text"/>
    <xsl:param name="from"/>
    <xsl:param name="to"/>
    <xsl:choose>
      <xsl:when test="contains($text,$from)">
        <xsl:value-of select="concat(substring-before($text,$from),$to)"/>
        <xsl:call-template name="str_replace">
          <xsl:with-param name="text" select="substring-after($text,$from)"/>
          <xsl:with-param name="from" select="$from"/>
          <xsl:with-param name="to"   select="$to"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$text"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template match="title">
    <xsl:call-template name="str_replace">
      <xsl:with-param name="text">
        <xsl:apply-templates/>
      </xsl:with-param>
	<xsl:with-param name="from" select="'_'"/>
	<xsl:with-param name="to" select="'\_'"/>
    </xsl:call-template>
  </xsl:template>

  <xsl:template match="title" name="title-strip" mode="strip">
    <xsl:value-of select="translate(.,'_','')"/>
  </xsl:template>


<!-- ========= obsolete_functions  and  obsolescent_functions ========== -->
<xsl:template match="obsolete_functions">
\begin{center}
\begin{longtable}{ll}
    <xsl:for-each select="//command[contains(title,'OBSOLETE')]">
{\verb~<xsl:value-of select="title"/>~} &amp;
      <xsl:apply-templates select="short_description"/>\\
   </xsl:for-each>
\end{longtable}
\end{center}
\noindent
</xsl:template>

<xsl:template match="obsolescent_functions">
\begin{center}
\begin{longtable}{ll}
    <xsl:for-each select="//command[contains(title,'OBSOLESCENT')]">
{\verb~<xsl:value-of select="title"/>~} &amp;
      <xsl:apply-templates select="short_description"/>\\
   </xsl:for-each>
\end{longtable}
\end{center}
\noindent
</xsl:template>


  <xsl:template match="commands_and_functions_for">  
   <xsl:variable name="my_type" select="@type"/>
\begin{center}
\begin{longtable}{ll}
   <xsl:for-each select="//command[type=$my_type]|//command[types/type=$my_type]|//command[syntax/type=$my_type]">
     <xsl:sort select="title"/>
{\verb~<xsl:value-of select="title"/>~} &amp;
      <xsl:apply-templates select="short_description"/>\\
   </xsl:for-each>
\end{longtable}
\end{center}

\noindent
</xsl:template>

  <xsl:template match="commands_and_functions_rtn">  
   <xsl:variable name="my_type" select="@type"/>
\begin{center}
\begin{longtable}{ll}
   <xsl:for-each select="//command[syntax/rtn=$my_type]">
     <xsl:sort select="title"/>
{\verb~<xsl:value-of select="title"/>~} &amp;
      <xsl:apply-templates select="short_description"/>\\
   </xsl:for-each>
\end{longtable}
\end{center}

\noindent
</xsl:template>

  <xsl:template match="coc_version">
    <xsl:value-of select="/help/version/cocoa_version"/>
  </xsl:template>

  <xsl:template match="i">\textit{<xsl:apply-templates/>}</xsl:template>
  <xsl:template match="b">\textbf{<xsl:apply-templates/>}</xsl:template>
  <xsl:template match="em">\textbf{<xsl:apply-templates/>}</xsl:template>

  <xsl:template match="quotes">&quot;<xsl:apply-templates/>&quot;</xsl:template>
  <xsl:template match="sup">^{<xsl:apply-templates/>}</xsl:template>
  <xsl:template match="tt" mode="verbatim">&quot;<xsl:apply-templates/>&quot;</xsl:template>  
  <xsl:template match="tt">``\verb&amp;<xsl:apply-templates/>&amp;''</xsl:template>

  <xsl:template match="formula">$<xsl:apply-templates/>$</xsl:template>

  <xsl:template match="less_eq">\le </xsl:template>
  <xsl:template match="times">\times </xsl:template>

  <xsl:template
  match="example">\begin{Verbatim}[label=example, rulecolor=\color{PineGreen}, frame=single]<xsl:apply-templates
  select="." mode="chooseverbatim"/>\end{Verbatim}
</xsl:template>

  <xsl:template match="par">\par </xsl:template>
  <xsl:template match="ie">\textit{i.e.}</xsl:template>
  <xsl:template match="eg">\textit{e.g.}</xsl:template>
  <xsl:template match="BOOL">\textsf{BOOL}</xsl:template>
  <xsl:template match="ERROR">\textsf{ERROR}</xsl:template>
  <xsl:template match="FUNCTION">\textsf{FUNCTION}</xsl:template>
  <xsl:template match="IDEAL">\textsf{IDEAL}</xsl:template>
  <xsl:template match="INT">\textsf{INT}</xsl:template>
  <xsl:template match="LIST">\textsf{LIST}</xsl:template>
  <xsl:template match="MAT">\textsf{MAT}</xsl:template>
  <xsl:template match="MODULE">\textsf{MODULE}</xsl:template>
  <xsl:template match="MODULEELEM">\textsf{MODULEELEM}</xsl:template>
  <xsl:template match="RAT">\textsf{RAT}</xsl:template>
  <xsl:template match="RECORD">\textsf{RECORD}</xsl:template>
  <xsl:template match="RING">\textsf{RING}</xsl:template>
  <xsl:template match="RINGELEM">\textsf{RINGELEM}</xsl:template>
  <xsl:template match="RINGHOM">\textsf{RINGHOM}</xsl:template>
  <xsl:template match="STRING">\textsf{STRING}</xsl:template>
  <xsl:template match="TAGGED">\textsf{TAGGED}</xsl:template>
  <xsl:template match="description/TYPE">\textsf{TYPE}</xsl:template>

  <xsl:template match="ref">``<xsl:apply-templates/>'' (\ref{<xsl:call-template name="title-strip" select="."/>} pg.\pageref{<xsl:call-template name="title-strip" select="."/>})</xsl:template>

  <xsl:template match="ttref">``\verb&amp;<xsl:apply-templates/>&amp;'' (\ref{<xsl:call-template name="title-strip" select="."/>} pg.\pageref{<xsl:call-template name="title-strip" select="."/>})</xsl:template>


  <xsl:template match="coc_date">
    <xsl:value-of select="/help/version/date"/>
  </xsl:template>

  <xsl:template match="verbatim">\begin{verbatim}<xsl:apply-templates
  select="." mode="chooseverbatim"/>\end{verbatim}</xsl:template>

<!--   <xsl:template match="list">\begin{itemize}<xsl:apply-templates -->
<!--   select="."/>\end{itemize}</xsl:template> -->

<!--   <xsl:template match="item">\item <xsl:apply-templates -->
<!--   select="."/></xsl:template> -->

  <xsl:template match="*" mode="chooseverbatim">
    <xsl:choose>
      <xsl:when test="verbatim='yes'">
        <xsl:apply-templates mode="verbatim"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:apply-templates/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

</xsl:stylesheet>

