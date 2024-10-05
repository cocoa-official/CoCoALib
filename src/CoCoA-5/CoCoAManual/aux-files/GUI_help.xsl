<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version="1.1"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:preserve-space elements="quotes"/>
<xsl:preserve-space elements="tt"/>

<!-- xslt document for creating the GUI-man directory for CoCoA -->

<xsl:output method="text"/>
<xsl:strip-space elements="*"/>

<!-- line breaks -->

<xsl:variable name="endl">
<xsl:text>
</xsl:text>
</xsl:variable>

<xsl:param name="CoCoAVersion" required="yes"/>

<xsl:template match="help">


<!-- ____ COMMANDS ______________________________________________ -->


  <xsl:for-each select="cocoa_commands/chapter_letter">
    <xsl:for-each select="command">
      <xsl:sort select="title"/>

      <xsl:variable name="file_cmd">
        <xsl:call-template name="create_href">
          <xsl:with-param name="key" select="title"/>
        </xsl:call-template>
      </xsl:variable>
      <xsl:document href="{$file_cmd}">

<xsl:text>
&lt;html&gt;
&lt;head&gt;
&lt;link rel="stylesheet" type="text/css" href="gui.css"&gt;
</xsl:text>

        <xsl:text>&lt;title&gt;</xsl:text>
        <xsl:apply-templates select="title"/>
        <xsl:text>&lt;/title&gt;</xsl:text>

<xsl:text>
&lt;/head&gt;

&lt;body bgcolor=#eeffff&gt;
&lt;div&gt;
</xsl:text>

        <xsl:text>&lt;a href=&quot;toc.html#</xsl:text>
        <xsl:value-of select="title"/>
        <xsl:text>&quot;&gt;up&lt;/a&gt;</xsl:text>
        <xsl:text> </xsl:text>

        <xsl:choose>
          <xsl:when test="position()!=1">
            <xsl:text>&lt;a href=&quot;</xsl:text>
            <xsl:call-template name="previous_href">
              <xsl:with-param name="key" select="title"/>
            </xsl:call-template>
            <xsl:text>&quot;&gt;previous&lt;/a&gt; </xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>&lt;span class="grayed"&gt;previous&lt;/span&gt;</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
        <xsl:text> </xsl:text>

        <xsl:choose>
          <xsl:when test="position()!=last()">
            <xsl:text>&lt;a href=&quot;</xsl:text>
              <xsl:call-template name="next_href">
              <xsl:with-param name="key" select="title"/>
            </xsl:call-template>
            <xsl:text>&quot;&gt;next&lt;/A&gt;</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>&lt;span class="grayed"&gt;next&lt;/span&gt;</xsl:text>
          </xsl:otherwise>
        </xsl:choose>

        <xsl:call-template name="h1">
          <xsl:with-param name="header">
	    <xsl:apply-templates select="title"/>
            <xsl:text>&lt;small&gt; &lt;i&gt;</xsl:text>
            <xsl:text> &amp;nbsp;&amp;nbsp; -- &amp;nbsp;&amp;nbsp; </xsl:text>
            <xsl:value-of select="$endl"/>
            <xsl:apply-templates select="short_description"/>
            <xsl:text>&lt;/i&gt; &lt;/small&gt;</xsl:text>
	  </xsl:with-param>
        </xsl:call-template>


<!--         <xsl:text>&lt;br&gt;</xsl:text> -->
<!--         <xsl:value-of select="$endl"/> -->
<!--         <xsl:apply-templates select="short_description"/> -->

	<xsl:apply-templates select="syntax"/>

        <xsl:call-template name="h2">
          <xsl:with-param name="header">Description</xsl:with-param>
        </xsl:call-template>
        <xsl:apply-templates select="description"/>

        <!-- new style "See also" -->
        <xsl:if test="seealso">
          <xsl:apply-templates select="seealso"/>
        </xsl:if>

        <!-- old style "See also" -->
<!--         <xsl:if test="see"> -->
<!--           <xsl:call-template name="seealso"/> -->
<!--         </xsl:if> -->

<xsl:text>
&lt;/div&gt;

&lt;/body&gt;
&lt;/html&gt;
</xsl:text>

      </xsl:document>
    </xsl:for-each>
  </xsl:for-each>

<!-- ____ CHAPTERS ______________________________________________ -->


  <xsl:for-each select="manual_parts/part">
    <xsl:variable name="part_num" select="position()"/>
    <xsl:for-each select="chapter">
       <xsl:variable name="chap_num" select="position()"/>
      <xsl:for-each select="section">
        <xsl:variable name="sect_num" select="position()"/>

        <xsl:variable name="file">
          <xsl:call-template name="create_href"><xsl:with-param name="key" select="title"/></xsl:call-template>
        </xsl:variable>

        <xsl:document href="{$file}">

<xsl:text>
&lt;html&gt;
&lt;head&gt;
&lt;link REL="stylesheet" TYPE="text/css" href="gui.css"&gt;
</xsl:text>

          <xsl:text>&lt;title&gt;</xsl:text>
          <xsl:apply-templates select="title"/>
          <xsl:text>&lt;/title&gt;</xsl:text>

<xsl:text>
&lt;/head&gt;

&lt;body bgcolor=#eeffff&gt;
&lt;div&gt;
</xsl:text>

          <xsl:text>&lt;a href=&quot;toc.html#</xsl:text>
          <xsl:value-of select="concat('toc.html#p',$part_num,'c',$chap_num)"/>
          <xsl:text>&quot;&gt;up&lt;/a&gt;</xsl:text>
          <xsl:text> </xsl:text>

          <xsl:choose>
            <xsl:when test="position()!=1">
              <xsl:text>&lt;a href=&quot;</xsl:text>
              <xsl:call-template name="previous_href">
                <xsl:with-param name="key" select="title"/>
              </xsl:call-template>
              <xsl:text>&quot;&gt;previous&lt;/a&gt;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>&lt;span class="grayed"&gt;previous&lt;/span&gt;</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
          <xsl:text> </xsl:text>

          <xsl:choose>
            <xsl:when test="position()!=last()">
              <xsl:text>&lt;a href=&quot;</xsl:text>
              <xsl:call-template name="next_href">
                <xsl:with-param name="key" select="title"/>
              </xsl:call-template>
              <xsl:text>&quot;&gt;next&lt;/a&gt;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>&lt;span class="grayed"&gt;next&lt;/span&gt;</xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <xsl:call-template name="h1">
            <xsl:with-param name="header">
              <xsl:value-of select="concat($part_num,'.',$chap_num,'.',$sect_num,' ')"/>
              <xsl:apply-templates select="title"/>
            </xsl:with-param>
          </xsl:call-template>

          <xsl:apply-templates select="description"/>

<xsl:text>
&lt;/div&gt;

&lt;/body&gt;
&lt;/html&gt;
</xsl:text>

        </xsl:document>
      </xsl:for-each>
    </xsl:for-each>
  </xsl:for-each>


<!-- ____ TOC ___________________________ -->


  <xsl:variable name="file_toc" select="'toc.html'"/>
  <xsl:document href="{$file_toc}">

<xsl:text>
&lt;HTML&gt;
&lt;HEAD&gt;
&lt;link REL="stylesheet" TYPE="text/css" href="gui.css"&gt;
</xsl:text>

<xsl:text>&lt;title&gt;</xsl:text>
<xsl:text>CoCoA: table of contents</xsl:text>
<xsl:text>&lt;/title&gt;</xsl:text>

<xsl:text>
&lt;/head&gt;

&lt;body bgcolor=#eeffff&gt;

&lt;a name=&quot;top&quot;&gt;&lt;/a&gt;

&lt;h1&gt;
</xsl:text>
  CoCoA-<xsl:value-of select="$CoCoAVersion"/> Manual
<xsl:text>
&lt;/h1&gt;
</xsl:text>

<xsl:call-template name="h2">
  <xsl:with-param name="header">
    <xsl:text>
&lt;a href=&quot;#alphalist&quot;&gt;Part 1: Alphabetical list of commands&lt;/a&gt;
&lt;br&gt;
&lt;a href=&quot;#2&quot;&gt;Part 2: CoCoA Language &amp; Tutorials&lt;/a&gt;
&lt;br&gt;
&lt;a href=&quot;#3&quot;&gt;Part 3: CoCoA Data Types&lt;/a&gt;
    </xsl:text>
  </xsl:with-param>
</xsl:call-template>


<!-- ____ alphabetical list of COMMANDS ___________________________ -->

<xsl:text>
&lt;a name=&quot;alphalist&quot;&gt;&lt;/a&gt;
&lt;div class=commands&gt;
</xsl:text>

    <xsl:call-template name="h1">
      <xsl:with-param name="header">
        <xsl:text>Part 1: Alphabetical list of commands</xsl:text>
      </xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="h2">
      <xsl:with-param name="header">
        <xsl:for-each select="/help/cocoa_commands/chapter_letter">
          <xsl:if test="not(string-length(title)>1)">
            <xsl:text>&lt;a href=&quot;#</xsl:text>
            <xsl:value-of select="title"/>
            <xsl:text>&quot;&gt;</xsl:text>
            <xsl:value-of select="title"/>
            <xsl:text>&lt;/a&gt; </xsl:text>
          </xsl:if>
        </xsl:for-each>
      </xsl:with-param>
    </xsl:call-template>

    <xsl:for-each select="cocoa_commands/chapter_letter">

      <xsl:text>&lt;a name=&quot;</xsl:text>
      <xsl:value-of select="title"/>
      <xsl:text>&quot;&gt;&lt;h3&gt;</xsl:text>
      <xsl:value-of select="title"/>
      <xsl:text>&lt;/h3&gt;&lt;/a&gt;</xsl:text>
      <xsl:value-of select="$endl"/>

      <xsl:text>&lt;table bgcolor=#ccffff width=100%&gt;&lt;tr&gt;&lt;td&gt;</xsl:text>

      <xsl:for-each select="command">
        <xsl:sort select="title"/>

        <xsl:value-of select="$endl"/>
        <xsl:text>&lt;a name=&quot;</xsl:text>
        <xsl:value-of select="title"/>
        <xsl:text>&quot;&gt;&lt;/a&gt;</xsl:text>
        <xsl:value-of select="$endl"/>

        <xsl:text>&lt;a href=&quot;</xsl:text>
        <xsl:call-template name="create_href">
          <xsl:with-param name="key" select="title"/>
        </xsl:call-template>
        <xsl:text>&quot;&gt;</xsl:text>

        <xsl:text>&lt;b&gt;</xsl:text>
        <xsl:apply-templates select="title"/>
        <xsl:text>&lt;/b&gt;&lt;/a&gt;</xsl:text>

        <xsl:text> -- </xsl:text>

        <xsl:apply-templates select="short_description"/>

        <xsl:value-of select="$endl"/>
        <xsl:text>&lt;br&gt;</xsl:text>

      </xsl:for-each>

<xsl:text>
&lt;/td&gt;&lt;/tr&gt;&lt;/table&gt;

&lt;br&gt;
&lt;i&gt;
&lt;a href=&quot;#top&quot;&gt;Back to the top&lt;/a&gt; ---
&lt;a href=&quot;#alphalist&quot;&gt;Back to the alphabet&lt;/a&gt;
&lt;/i&gt;
&lt;br&gt;
</xsl:text>

    </xsl:for-each>

<!-- ____ table of contents CHAPTERS ____________________________ -->

<!-- enable/disable toc about chapters -->
  <xsl:if test="true()">
    <xsl:for-each select="/help/manual_parts/part">
      <xsl:variable name="part_num" select="position()+1"/>

      <xsl:text>&lt;div class=p</xsl:text>
      <xsl:value-of select="$part_num"/>
      <xsl:text>&gt;</xsl:text>
      <xsl:value-of select="$endl"/>
      <xsl:text>&lt;a name=&quot;</xsl:text>
      <xsl:value-of select="$part_num"/>
      <xsl:text>&quot;&gt;&lt;/a&gt;</xsl:text>

      <xsl:call-template name="h1">
        <xsl:with-param name="header">
          <xsl:text>Part </xsl:text><xsl:value-of select="$part_num"/>
          <xsl:text>:  </xsl:text><xsl:value-of select="title"/>
        </xsl:with-param>
      </xsl:call-template>

      <xsl:for-each select="chapter">
        <xsl:variable name="chap_num" select="position()"/>

        <xsl:text>&lt;A NAME=&quot;</xsl:text>
        <xsl:value-of select="concat('p',$part_num,'c',$chap_num)"/>
        <xsl:text>&quot;&gt;&lt;/A&gt;</xsl:text>
        <xsl:text>&lt;A NAME=&quot;</xsl:text>
        <xsl:apply-templates select="title"/>
        <xsl:text>&quot;&gt;&lt;/A&gt;</xsl:text>

        <xsl:call-template name="h3">
          <xsl:with-param name="header">
            <xsl:text>Part </xsl:text><xsl:value-of select="$part_num"/>
            <xsl:text> - </xsl:text>
            <xsl:text>Chapter </xsl:text><xsl:value-of select="$chap_num"/>
            <xsl:text> - </xsl:text>
            <xsl:apply-templates select="title"/>
          </xsl:with-param>
        </xsl:call-template>

        <xsl:value-of select="$endl"/>
        <xsl:text>&lt;table bgcolor=#ccffff width=100%&gt;&lt;tr&gt;&lt;td&gt;</xsl:text>

        <xsl:for-each select="section">
          <xsl:variable name="sect_num" select="position()"/>

          <xsl:value-of select="$endl"/>
          <xsl:text>&lt;A NAME=&quot;</xsl:text>
          <xsl:value-of select="title"/>
          <xsl:text>&quot;&gt;&lt;/A&gt;</xsl:text>
          <xsl:value-of select="$endl"/>

          <xsl:value-of select="$part_num"/><xsl:text>.</xsl:text>
          <xsl:value-of select="$chap_num"/><xsl:text>.</xsl:text>
          <xsl:value-of select="$sect_num"/><xsl:text>.</xsl:text>

          <xsl:text> &lt;A HREF=&quot;</xsl:text>
          <xsl:call-template name="create_href">
            <xsl:with-param name="key" select="title"/>
          </xsl:call-template>
          <xsl:text>&quot;&gt;</xsl:text>
          <xsl:apply-templates select="title"/>
          <xsl:text>&lt;/A&gt;</xsl:text>

          <xsl:value-of select="$endl"/>
          <xsl:text>&lt;BR&gt;</xsl:text>

        </xsl:for-each>

<xsl:text>
&lt;/td&gt;&lt;/tr&gt;&lt;/table&gt;
&lt;I&gt;&lt;A HREF=&quot;#top&quot;&gt;Back to the top&lt;/A&gt;&lt;/I&gt;
&lt;BR&gt;
</xsl:text>

      </xsl:for-each>

<xsl:text>
&lt;/DIV&gt;
&lt;BR&gt;
</xsl:text>

    </xsl:for-each>
</xsl:if>
<!-- enable/disable create toc about chapters -->


<xsl:text>
&lt;/DIV&gt;
&lt;/BODY&gt;
&lt;/HTML&gt;
</xsl:text>

  </xsl:document>
</xsl:template>


<!-- ____ TEMPLATES ______________________________________________ -->

<xsl:template match="seealso" name="seealso">
          <xsl:call-template name="h2">
            <xsl:with-param name="header">See Also</xsl:with-param>
          </xsl:call-template>

          <xsl:text>&lt;ul&gt;</xsl:text>
          <xsl:value-of select="$endl"/>

          <xsl:for-each select="see">
            <xsl:text>&lt;li&gt;&lt;a href=&quot;</xsl:text>
            <xsl:call-template name="create_href">
              <xsl:with-param name="key" select="."/>
            </xsl:call-template>
            <xsl:text>&quot;&gt;</xsl:text>
            <xsl:value-of select="."/>
            <xsl:text>&lt;/A&gt;</xsl:text>
            <xsl:value-of select="$endl"/>
          </xsl:for-each>

          <xsl:text>&lt;/ul&gt;</xsl:text>
</xsl:template>

<!-- ========= obsolete_functions  and  obsolescent_functions ========== -->
<xsl:template match="obsolete_functions">
<xsl:text>
&lt;CENTER&gt;
&lt;TABLE BORDER=1&gt;
</xsl:text>

<xsl:for-each select="//command[contains(title,'OBSOLETE')]">
<xsl:text>
  &lt;TR bgcolor=#ffffff&gt;&lt;TD&gt;
</xsl:text>
    <xsl:text>&lt;A HREF=&quot;</xsl:text>
    <xsl:call-template name="create_href">
      <xsl:with-param name="key" select="title"/>
    </xsl:call-template>
    <xsl:text>&quot;&gt;</xsl:text>
    <xsl:value-of select="title"/>
    <xsl:text>&lt;/A&gt;</xsl:text>
<xsl:text>
  &lt;/TD&gt;&lt;TD&gt;
</xsl:text>
<xsl:apply-templates select="short_description"/>
<xsl:text>
  &lt;/TD&gt;&lt;/TR&gt;
</xsl:text>

</xsl:for-each>

<xsl:text>
&lt;/TABLE&gt;
&lt;/CENTER&gt;
&lt;BR&gt;&lt;BR&gt;
</xsl:text>
</xsl:template>

<!-- ===== -->
<xsl:template match="obsolescent_functions">
<xsl:text>
&lt;CENTER&gt;
&lt;TABLE BORDER=1&gt;
</xsl:text>

<xsl:for-each select="//command[contains(title,'OBSOLESCENT')]">
<xsl:text>
  &lt;TR bgcolor=#ffffff&gt;&lt;TD&gt;
</xsl:text>
    <xsl:text>&lt;A HREF=&quot;</xsl:text>
    <xsl:call-template name="create_href">
      <xsl:with-param name="key" select="title"/>
    </xsl:call-template>
    <xsl:text>&quot;&gt;</xsl:text>
    <xsl:value-of select="title"/>
    <xsl:text>&lt;/A&gt;</xsl:text>
<xsl:text>
  &lt;/TD&gt;&lt;TD&gt;
</xsl:text>
<xsl:apply-templates select="short_description"/>
<xsl:text>
  &lt;/TD&gt;&lt;/TR&gt;
</xsl:text>
</xsl:for-each>

<xsl:text>
&lt;/TABLE&gt;
&lt;/CENTER&gt;
&lt;BR&gt;&lt;BR&gt;
</xsl:text>
</xsl:template>

<!-- ========= commands_and_functions_for ========== -->
<xsl:template match="commands_and_functions_for">
<xsl:text>
&lt;CENTER&gt;
&lt;TABLE BORDER=1&gt;
</xsl:text>

<xsl:variable name="my_type" select="@type"/>
<xsl:for-each select="//command[type=$my_type]|//command[types/type=$my_type]|//command[syntax/type=$my_type]">
  <xsl:sort select="title"/>

<xsl:text>
  &lt;TR bgcolor=#ffffff&gt;&lt;TD&gt;
</xsl:text>
    <xsl:text>&lt;A HREF=&quot;</xsl:text>
    <xsl:call-template name="create_href">
      <xsl:with-param name="key" select="title"/>
    </xsl:call-template>
    <xsl:text>&quot;&gt;</xsl:text>
    <xsl:value-of select="title"/>
    <xsl:text>&lt;/A&gt;</xsl:text>
<xsl:text>
  &lt;/TD&gt;&lt;TD&gt;
</xsl:text>
<xsl:apply-templates select="short_description"/>
<xsl:text>
  &lt;/TD&gt;&lt;/TR&gt;
</xsl:text>

</xsl:for-each>

<xsl:text>
&lt;/TABLE&gt;
&lt;/CENTER&gt;
&lt;BR&gt;&lt;BR&gt;
</xsl:text>
</xsl:template>

<!-- ========= commands_and_functions_rtn ========== -->
<xsl:template match="commands_and_functions_rtn">
<xsl:text>
&lt;CENTER&gt;
&lt;TABLE BORDER=1&gt;
</xsl:text>

<xsl:variable name="my_type" select="@type"/>
<xsl:for-each select="//command[syntax/rtn=$my_type]">
  <xsl:sort select="title"/>

<xsl:text>
  &lt;TR bgcolor=#ffffff&gt;&lt;TD&gt;
</xsl:text>
    <xsl:text>&lt;A HREF=&quot;</xsl:text>
    <xsl:call-template name="create_href">
      <xsl:with-param name="key" select="title"/>
    </xsl:call-template>
    <xsl:text>&quot;&gt;</xsl:text>
    <xsl:value-of select="title"/>
    <xsl:text>&lt;/A&gt;</xsl:text>
<xsl:text>
  &lt;/TD&gt;&lt;TD&gt;
</xsl:text>
<xsl:apply-templates select="short_description"/>
<xsl:text>
  &lt;/TD&gt;&lt;/TR&gt;
</xsl:text>

</xsl:for-each>

<xsl:text>
&lt;/TABLE&gt;
&lt;/CENTER&gt;
&lt;BR&gt;&lt;BR&gt;
</xsl:text>
</xsl:template>
<!-- ========= end == commands_and_functions_rtn ========== -->

<xsl:template match="em">&lt;b&gt;<xsl:value-of select="."/>&lt;/b&gt;</xsl:template>
<xsl:template match="b">&lt;b&gt;<xsl:value-of select="."/>&lt;/b&gt;</xsl:template>
<xsl:template match="i">&lt;i&gt;<xsl:value-of select="."/>&lt;/i&gt;</xsl:template>

<xsl:template match="quotes">&quot;<xsl:apply-templates/>&quot;</xsl:template>
<xsl:template match="sup">&lt;sup&gt;<xsl:apply-templates/>&lt;/sup&gt;</xsl:template>
<xsl:template match="formula">&lt;font color=#0000aa&gt;<xsl:apply-templates/>&lt;/font&gt;</xsl:template>

<xsl:template match="tt">&lt;tt&gt;&lt;b style="background:#bbeeee"&gt;<xsl:apply-templates/>&lt;/b&gt;&lt;/tt&gt;</xsl:template>

<!-- <xsl:template match="backslash">\<xsl:apply-templates/></xsl:template> -->
<xsl:template match="less_eq">&amp;lt;=</xsl:template>
<xsl:template match="times"> x </xsl:template>

<!--   <xsl:template match="ref">&lt;A HREF=&quot;toc.html#<xsl:apply-templates/>&quot;&gt;&quot;<xsl:apply-templates/>&quot;&lt;/A&gt;</xsl:template> -->

<xsl:template match="ref">
  <xsl:text>&lt;A HREF=&quot;</xsl:text>
  <xsl:call-template name="create_href">
    <xsl:with-param name="key" select="."/>
  </xsl:call-template>
  <xsl:text>&quot;&gt;</xsl:text>
  <xsl:value-of select="."/>
  <xsl:text>&lt;/A&gt;</xsl:text>
</xsl:template>


<xsl:template match="ttref">
  <xsl:text>&lt;A HREF=&quot;</xsl:text>
  <xsl:call-template name="create_href">
    <xsl:with-param name="key" select="."/>
  </xsl:call-template>
  <xsl:text>&quot;&gt;</xsl:text>
    &lt;tt&gt;<xsl:value-of select="."/>&lt;/tt&gt;
  <xsl:text>&lt;/A&gt;</xsl:text>
</xsl:template>


<xsl:template match="par">&lt;br&gt;&lt;br&gt;</xsl:template>
<xsl:template match="ie">i.e.</xsl:template>
<xsl:template match="eg">e.g.</xsl:template>
<xsl:template match="BOOL">BOOL</xsl:template>
<xsl:template match="ERROR">ERROR</xsl:template>
<xsl:template match="FUNCTION">FUNCTION</xsl:template>
<xsl:template match="IDEAL">IDEAL</xsl:template>
<xsl:template match="INT">INT</xsl:template>
<xsl:template match="LIST">LIST</xsl:template>
<xsl:template match="MAT">MAT</xsl:template>
<xsl:template match="MODULE">MODULE</xsl:template>
<xsl:template match="MODULEELEM">MODULEELEM</xsl:template>
<xsl:template match="RAT">RAT</xsl:template>
<xsl:template match="RECORD">RECORD</xsl:template>
<xsl:template match="RING">RING</xsl:template>
<xsl:template match="RINGELEM">RINGELEM</xsl:template>
<xsl:template match="RINGHOM">RINGHOM</xsl:template>
<xsl:template match="STRING">STRING</xsl:template>
<xsl:template match="TAGGED">TAGGED</xsl:template>
<xsl:template match="description/TYPE">TYPE</xsl:template>

<xsl:template name="header_and_verbatim">
  <xsl:param name="header"/>
  <xsl:text>&lt;br&gt;</xsl:text>
  <xsl:call-template name="h3">
    <xsl:with-param name="header"><xsl:value-of select="$header"/></xsl:with-param>
  </xsl:call-template>
  <xsl:text>&lt;table bgcolor=#ccffff width=100%&gt;</xsl:text>
  <xsl:text>&lt;tr&gt;&lt;td&gt;</xsl:text>
  <xsl:text>&lt;pre&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/pre&gt;</xsl:text>
  <xsl:text>&lt;/td&gt;&lt;/tr&gt;</xsl:text><!-- /td/tr -->
  <xsl:text>&lt;/table&gt;</xsl:text><!-- /table -->
  <xsl:text>
</xsl:text><!-- newline -->
</xsl:template>


<xsl:template match="example">
  <xsl:call-template name="header_and_verbatim">
    <xsl:with-param name="header">Example</xsl:with-param>
  </xsl:call-template>
</xsl:template>


<xsl:template match="syntax">
  <xsl:call-template name="header_and_verbatim">
    <xsl:with-param name="header">Syntax</xsl:with-param>
  </xsl:call-template>
</xsl:template>


<xsl:template match="verbatim">
  &lt;pre&gt;<xsl:apply-templates/>&lt;/pre&gt;
</xsl:template>


<xsl:template name="hN">
  <xsl:param name="header" select="A"/>
  <xsl:param name="font_size" select="+1"/>
  <xsl:value-of select="$endl"/>
  <xsl:text>
</xsl:text><!-- newline -->
  <xsl:text>&lt;br&gt; </xsl:text>
  <xsl:text>
</xsl:text><!-- newline -->
  <xsl:text>&lt;!-- ========================= --&gt;</xsl:text>
  <xsl:text>
</xsl:text><!-- newline -->
  <xsl:text>&lt;table bgcolor=#00dddd width=100%&gt;</xsl:text>
  <xsl:text>&lt;tr&gt;&lt;td&gt;</xsl:text><!-- tr td -->
  <xsl:value-of select="$endl"/>
  <xsl:text>  &lt;font size=</xsl:text>
  <xsl:value-of select="$font_size"/>
  <xsl:text>&gt;&lt;b&gt;</xsl:text><!-- font..b -->
  <xsl:value-of select="$header"/>
  <xsl:text>&lt;/b&gt;&lt;/font&gt;</xsl:text><!-- /font../b -->
  <xsl:value-of select="$endl"/>
  <xsl:text>&lt;/td&gt;&lt;/tr&gt;&lt;/table&gt;</xsl:text><!-- /tr/td/table-->
  <xsl:text>
</xsl:text><!-- newline -->
</xsl:template>


<xsl:template name="h1">
  <xsl:param name="header" select="A"/>
  <xsl:call-template name="hN">
    <xsl:with-param name="header"><xsl:value-of select="$header"/></xsl:with-param>
    <xsl:with-param name="font_size">+3</xsl:with-param>
  </xsl:call-template>
</xsl:template>


<xsl:template name="h2">
  <xsl:param name="header" select="A"/>
  <xsl:call-template name="hN">
    <xsl:with-param name="header"><xsl:value-of select="$header"/></xsl:with-param>
    <xsl:with-param name="font_size">+2</xsl:with-param>
  </xsl:call-template>
</xsl:template>


<xsl:template name="h3">
  <xsl:param name="header" select="A"/>
  <xsl:call-template name="hN">
    <xsl:with-param name="header"><xsl:value-of select="$header"/></xsl:with-param>
    <xsl:with-param name="font_size">+1</xsl:with-param>
  </xsl:call-template>
</xsl:template>


<!-- ____ direct hyperlink references ___________________________ -->

<!-- since we can't build associative arrays we'll have to collect key/value pairs for the links -->

<xsl:variable name="command_hrefs">
  <xsl:for-each select="//cocoa_commands/chapter_letter">
    <xsl:variable name="initial" select="title"/>
  <xsl:for-each select="command">
    <xsl:sort select="title"/>
    <link>
      <key><xsl:value-of select="title"/></key>
      <value>
        <xsl:choose>
          <xsl:when test="$initial != 'Symbol' and function-available('translate')">
            <xsl:value-of select="concat('cmd',translate(title,',.:/()- ',''),'.html')"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="concat('cmd',$initial,position(),'.html')"/>
          </xsl:otherwise>
        </xsl:choose>
      </value>
    </link>
  </xsl:for-each>
  </xsl:for-each>
</xsl:variable>


<xsl:variable name="parts_hrefs">
  <xsl:for-each select="//manual_parts/part">
    <xsl:variable name="part_num" select="position()"/>
    <xsl:for-each select="chapter">
      <xsl:variable name="chap_num" select="position()"/>
      <xsl:for-each select="section">
        <xsl:variable name="sect_num" select="position()"/>
        <link>
          <key><xsl:value-of select="title"/></key>
          <value>
            <xsl:choose>
              <xsl:when test="function-available('translate')">
                <xsl:value-of select="concat('part',translate(title,',.:/()- ',''),'.html')"/>
              </xsl:when>
              <xsl:otherwise>
                <xsl:value-of select="concat('p',$part_num,'c',$chap_num, 's',$sect_num,'.html')"/>
              </xsl:otherwise>
            </xsl:choose>
          </value>
        </link>
      </xsl:for-each>
    </xsl:for-each>
  </xsl:for-each>
</xsl:variable>


<xsl:template name="create_href">
  <xsl:param name="key"/> <!-- it is assumed that the key is unique... -->
  <xsl:choose>
    <xsl:when test="$command_hrefs/*[key=$key]">
      <xsl:for-each select="$command_hrefs/*[key=$key]">
        <xsl:value-of select="value"/>
      </xsl:for-each>
    </xsl:when>
    <xsl:when test="$parts_hrefs/*[key=$key]">
      <xsl:for-each select="$parts_hrefs/*[key=$key]">
        <xsl:value-of select="value"/>
      </xsl:for-each>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>toc.html#</xsl:text>
      <xsl:value-of select="$key"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>


<xsl:template name="previous_href">
  <xsl:param name="key"/>
  <xsl:choose>
    <xsl:when test="$command_hrefs/*[key=$key]">
      <xsl:for-each select="$command_hrefs/*[key=$key]">
        <xsl:value-of select="preceding-sibling::*[1]/value"/>
      </xsl:for-each>
    </xsl:when>
    <xsl:when test="$parts_hrefs/*[key=$key]">
      <xsl:for-each select="$parts_hrefs/*[key=$key]">
        <xsl:value-of select="preceding-sibling::*[1]/value"/>
      </xsl:for-each>
    </xsl:when>
  </xsl:choose>
</xsl:template>


<xsl:template name="next_href">
  <xsl:param name="key"/>
  <xsl:choose>
    <xsl:when test="$command_hrefs/*[key=$key]">
      <xsl:for-each select="$command_hrefs/*[key=$key]">
        <xsl:value-of select="following-sibling::*[1]/value"/>
      </xsl:for-each>
    </xsl:when>
    <xsl:when test="$parts_hrefs/*[key=$key]">
      <xsl:for-each select="$parts_hrefs/*[key=$key]">
        <xsl:value-of select="following-sibling::*[1]/value"/>
      </xsl:for-each>
    </xsl:when>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
