﻿<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en" >
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
<meta name="robots" content="index,follow" />
<meta name="verify-v1" content="6esrL4Jqlv4/2lbCRFu0SeYW5zrKD7E9U+gmGlKNtro=" />
<META NAME="keywords" CONTENT="blast,basic local alignment search tool,sequence similarity,sequence alignment,sequence homology,bioinformatics" />
<META NAME="description" CONTENT="The Basic Local Alignment Search Tool (BLAST) finds regions of local similarity between sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance of matches. BLAST can be used to infer functional and evolutionary relationships between sequences as well as help identify members of gene families." />
<meta name="ncbitoggler" content="animation:'none'"/>
<meta name="referrer" content="origin-when-cross-origin" />
<meta name="ncbi_app" content="blast" />
<meta name="ncbi_pdid" content="blasthome" />

<meta name="ncbi_db" content="" />
<meta name="ncbi_program" content="" />
<meta name="ncbi_algorithm" content="" />

<meta name="ncbi_stat" content="false" />
<meta name="ncbi_sessionid" content="F4144FF8278C60E3_18068SID" />
<meta name="ncbi_phid" content="50C91C232B09E8610000000000000001" />
<meta name="ncbi_pcid" content="" />
<script type="text/javascript"> var ncbi_startTime = new Date(); </script>
<title>BLAST: Basic Local Alignment Search Tool</title>
<script type="text/javascript" src="/core/jig/1.15.2/js/jig.min.js             "></script>
<script type="text/javascript">    jQuery.getScript("/core/alerts/alerts.js", function() {
        galert(['div#header', 'body > *:nth-child(1)'])
    });</script>
<meta http-equiv="Pragma" content="no-cache">
<link rel="stylesheet" type="text/css" href="css/uswds.min.css" media="screen" />
<link rel="stylesheet"  type="text/css" href="https://www.ncbi.nlm.nih.gov/style-guide/static/nwds/css/nwds.css"/>
<link rel="stylesheet" href="css/headerNew.css?v=1"/>
<link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.5.0/css/all.css" crossorigin="anonymous"> <!-- Font Awesome icons -->
<link rel="stylesheet" type="text/css" href="css/footerNew.css?v=1" media="screen" />



<link rel="stylesheet" type="text/css" href="css/home.css" media="screen" />
<link rel="stylesheet" type="text/css" href="css/print.css" media="print" />
<!--[if lte IE 6]>
<link rel="stylesheet" type="text/css" href="css/ie6_or_less.css" />
<![endif]-->
<script type="text/javascript" src="js/utils.js"></script>
<script type="text/javascript" src="js/blast.js"></script>
<script type="text/javascript" src="js/remote_data_provider.js"></script>
<script type="text/javascript" src="js/home.js"></script>
<script type="text/javascript" src="js/blastSurvey.js"></script>
<meta name="google-site-verification" content="DNp69-s7ILfn5FWNLN_xE4zVhpY0dDtLrLyOB_WOzpA" />
</head>
<body id="type-d">
<!--<div id="browsers_ajax"></div>-->
<script>var useOfficialGovtHeader = true;</script>
<section class="usa-banner">
  <div class="usa-accordion">
    <header class="usa-banner-header">
      <div class="usa-grid usa-banner-inner">
        <img src="https://www.ncbi.nlm.nih.gov/coreutils/uswds/img/favicons/favicon-57.png" alt="U.S. flag">
        <p>An official website of the United States government</p>
        <button class="usa-accordion-button usa-banner-button" aria-expanded="false" aria-controls="gov-banner-top">
          <span class="usa-banner-button-text">Here's how you know</span>
        </button>
      </div>
    </header>
    <div class="usa-banner-content usa-grid usa-accordion-content" id="gov-banner-top" aria-hidden="true">
      <div class="usa-banner-guidance-gov usa-width-one-half">
        <img class="usa-banner-icon usa-media_block-img" src="https://www.ncbi.nlm.nih.gov/coreutils/uswds/img/icon-dot-gov.svg" alt="Dot gov">
        <div class="usa-media_block-body">
          <p>
            <strong>The .gov means it’s official.</strong>
            <br>
            Federal government websites often end in .gov or .mil. Before
            sharing sensitive information, make sure you’re on a federal
            government site.
          </p>
        </div>
      </div>
      <div class="usa-banner-guidance-ssl usa-width-one-half">
        <img class="usa-banner-icon usa-media_block-img" src="https://www.ncbi.nlm.nih.gov/coreutils/uswds/img/icon-https.svg" alt="Https">
        <div class="usa-media_block-body">
          <p>
            <strong>The site is secure.</strong>
            <br>
            The <strong>https://</strong> ensures that you are connecting to the
            official website and that any information you provide is encrypted
            and transmitted securely.
          </p>
        </div>
      </div>
    </div>
  </div>
</section>
    	
<header class="ncbi-header" role="banner" data-section="Header">
<a class="usa-skipnav" href="#mainCont">Skip to main page content</a>
<div class="usa-grid">
    <div class="usa-width-one-whole">
        <div class="ncbi-header__logo">
                <a href="https://www.ncbi.nlm.nih.gov/" class="logo" aria-label="NCBI Logo" data-ga-action="click_image" data-ga-label="NIH NLM Logo">
                  <img src="https://www.ncbi.nlm.nih.gov/coreutils/nwds/img/logos/AgencyLogo.svg" alt="NIH NLM Logo">
                </a>
            </div>

        <div class="ncbi-header__account">
            <a id="account_login" href="https://www.ncbi.nlm.nih.gov/account/?back_url=https%3A%2F%2Fblast%2Encbi%2Enlm%2Enih%2Egov%2FBlast%2Ecgi%3FBLAST%5FPROGRAMS%3D%26PAGE%5FTYPE%3DBlastHome" class="usa-button header-button">Log in</a>
            <button id="account_info" class="header-button" aria-controls="account_popup">
                <span class="fa fa-user" aria-hidden="true"></span>
                <span class="username desktop-only" aria-hidden="true" id="uname_short"></span>
                <span class="sr-only">Show account info</span>
            </button>
        </div>

        <div class="ncbi-popup-anchor">
            <div class="ncbi-popup account-popup" id="account_popup" aria-hidden="true" role="dialog" aria-labelledby="account-popup-header">
                <div class="ncbi-popup-head">
                    <button class="ncbi-close-button"><span class="fa fa-window-close"></span><span class="usa-sr-only">Close</span></button>
                    <h4>Account</h4>
                </div>
                <div class="account-user-info">
                    Logged in as:<br>
                    <b><span class="username" id="uname_long">username</span></b>
                </div>
                <div class="account-links">
                    <ul class="usa-unstyled-list">
                        <li><a id="account_myncbi" href="https://www.ncbi.nlm.nih.gov/myncbi/">Dashboard</a> <span class="ncbi-text-small-light">(My NCBI)</span></li>
                        <li><a id="account_pubs" href="https://www.ncbi.nlm.nih.gov/myncbi/collections/bibliography/">Publications</a> <span class="ncbi-text-small-light">(My Bibliography)</span></li>
                        <li><a id="account_settings" href="https://www.ncbi.nlm.nih.gov/account/settings/">Account settings</a></li>
                        <li><a id="account_logout" href="https://www.ncbi.nlm.nih.gov/account/signout/?back_url=https%3A%2F%2Fblast%2Encbi%2Enlm%2Enih%2Egov%2FBlast%2Ecgi%3FBLAST%5FPROGRAMS%3D%26PAGE%5FTYPE%3DBlastHome">Log out</a></li>
                    </ul>
                </div>
            </div>
        </div>

    </div>
</div>
</header>
<div role="navigation" aria-label="access keys">
<a id="nws_header_accesskey_0" href="https://www.ncbi.nlm.nih.gov/guide/browsers/#ncbi_accesskeys" class="usa-sr-only" accesskey="0" tabindex="-1">Access keys</a>
<a id="nws_header_accesskey_1" href="https://www.ncbi.nlm.nih.gov" class="usa-sr-only" accesskey="1" tabindex="-1">NCBI Homepage</a>
<a id="nws_header_accesskey_2" href="/myncbi/" class="set-base-url usa-sr-only" accesskey="2" tabindex="-1">MyNCBI Homepage</a>
<a id="nws_header_accesskey_3" href="#maincontent" class="usa-sr-only" accesskey="3" tabindex="-1">Main Content</a>
<a id="nws_header_accesskey_4" href="#" class="usa-sr-only" accesskey="4" tabindex="-1">Main Navigation</a>
</div>
<nav class="ncbi-topnav" id="navcontent">
    <div class="usa-grid">
        <a class="ncbi-topnav-root" href="Blast.cgi">BLAST <sup>&reg;</sup></a> <span id="brc"></span>
        <ul class="rf ncbi-topnav-list" id="topnav-list">
            <li class="first active"><a href="Blast.cgi?CMD=Web&amp;PAGE_TYPE=BlastHome" title="BLAST Home">Home</a></li>
            <li class="recent "><a href="Blast.cgi?CMD=GetSaved&amp;RECENT_RESULTS=on" title="Unexpired BLAST jobs">Recent Results</a></li>                
            <li class="saved "><a href="Blast.cgi?CMD=GetSaved" title="Saved sets of BLAST search parameters">Saved Strategies</a></li>
            <li  class= "last documentation "> <a href="Blast.cgi?CMD=Web&amp;PAGE_TYPE=BlastDocs" title="BLAST documentation">Help</a></li>                            
        </ul>
    </div>
</nav>


<div id="wrap">    
    <div id="content-wrap">
        <div id="content">
            <div id="top">
            <div id="homeDescr" class="clearfix">
            <h1 class="lf">Basic Local Alignment Search Tool</h1>	
            <div  id="blast_desc"><span><span>BLAST</span> finds regions of similarity between biological sequences. The program compares nucleotide or
            protein	sequences to sequence databases and calculates the statistical
            significance.
            <a id="moreHelp" class="helplink  ui-ncbitoggler" data-jig="ncbitoggler" data-jigconfig="isIcon:false" data-ncbitoggler-toggles="hpHelp" title="help" href="#">Learn more</a>
            </span>
            
            <div class="ui-helper-reset" aria-live="assertive" >
<div class="helpbox ui-ncbitoggler-slave" id="hpHelp">
The Basic Local Alignment Search Tool (BLAST) finds regions of	local
similarity between sequences. The program compares nucleotide or
protein	sequences to sequence databases and calculates the statistical
significance of matches. BLAST can be used to infer functional and
evolutionary relationships between sequences as well as help identify
members of gene families.
</div>
</div><!--ARIA-->
            </div>            
            <div class="newsbox">
            <!--<img alt="News" id="imNews" class="left norm_height" src="images/news-titlebar.png" />-->
            <div class="featurebox">
<h3>News</h3>
<dl>
<dt>BLAST+ 2.13.0 is here!</dt>

<dd>
<br><div><div class="ExternalClassA582325AABFF48CA9B186D881BF5CE2C"><p>Starting with this release, we are including the blastn_vdb and tblastn_vdb executables in the BLAST+ distribution.</p>

<br>
</div>Thu, 17 Mar 2022 12:00:00 EST</dd>
</dl>
<p><a class="morelink" href="Blast.cgi?CMD=Web&amp;PAGE_TYPE=BlastNews">More BLAST news...</a></p>
</div>           
            </div>
            

            </div>

            <div id="chooseprog" class="section clearfix">
            <h2>Web BLAST</h2>
            <a href="Blast.cgi?PROGRAM=blastn&amp;PAGE_TYPE=BlastSearch&amp;LINK_LOC=blasthome" class="left spFirst" id="homeBlastn" title="Nucleotide BLAST"><img alt="Nucleotide query->Nucleotide database" src="images/nucleutide-blast-cover.png"/></a>
            <div id="transl" class="left">
            <a href="Blast.cgi?PROGRAM=blastx&amp;PAGE_TYPE=BlastSearch&amp;LINK_LOC=blasthome" id="homeBlastx" title="blastx"><img class="" alt="Translated Nucleotide query->Protein database" src="images/blastx-arrow.png"/></a>
            <a href="Blast.cgi?PROGRAM=tblastn&amp;PAGE_TYPE=BlastSearch&amp;LINK_LOC=blasthome" id="homeTblastn" title="tblastn"><img class="" alt="Protein query->Translated Nucleotide database" src="images/tblastn-arrow.png"/></a>
            </div>
            <a href="Blast.cgi?PROGRAM=blastp&amp;PAGE_TYPE=BlastSearch&amp;LINK_LOC=blasthome" class="left" id="homeBlastp" title="Protein BLAST"><img alt="Protein query->Protein database" src="images/protein-blast-cover.png"/></a>
            </div>
            <div id="gn_search">
<form id="fgs" action="Blast.cgi" method="post">
<h3><label for="qorganism" class="bl" id="srchEuk">BLAST Genomes</label></h3>
<input name="TAXID" size="55"  type="text" topoption="" id="qorganism" value="" data-jigconfig="dictionary:'bdb_euk_ref_viruses_rep_genomes_sg',isCrossDomain:false" autocomplete="off" data-jig="ncbiautocomplete" class="reset bl" suggesthint="Enter organism common name, scientific name, or tax id"/>
<input id="srchDB" type="button" value="Search"/>
<input name="PAGE_TYPE" type="hidden" value="BlastSearch" />
<input name="SEARCH_INIT" type="hidden" value="ReprGenomeDBSearch" />
<input name="LINK_LOC" type="hidden" value="blasthomeDBSearch" />
<div id="dbStatInfo"></div>
<p class="help hidden" id="sgHelp">
<span id="suggestPrompt">
Enter organism common name, scientific name, or tax id.
</span>
</p>
</form>
<ul>
<li><a href="Blast.cgi?PAGE_TYPE=BlastSearch&amp;BLAST_SPEC=OGP__9606__9558&amp;LINK_LOC=blasthome" title="BLAST human genome">Human</a></li>
<li><a href="Blast.cgi?PAGE_TYPE=BlastSearch&amp;BLAST_SPEC=OGP__10090__9559&amp;LINK_LOC=blasthome" title="BLAST mouse genome">Mouse</a></li>
<li><a href="Blast.cgi?PAGE_TYPE=BlastSearch&amp;BLAST_SPEC=OGP__10116__10621&amp;LINK_LOC=blasthome" title="BLAST rat genome">Rat</a></li>
<li><a href="Blast.cgi?PAGE_TYPE=BlastSearch&amp;BLAST_SPEC=MicrobialGenomes" title="BLAST microbial genomes">Microbes</a></li>
</ul>
</div>
            </div>
            <div class="perksblast section">
<h2>Standalone and API BLAST</h2>
<ul>
<li>
<a class="prksLink" href="Blast.cgi?CMD=Web&amp;PAGE_TYPE=BlastDocs&amp;DOC_TYPE=Download" title="Download Blast">
    <div id="dwnld-bg-img" class="bg-img">        
         <div class="prksTitle">Download BLAST</div>
         <div class="prksDescr">Get BLAST databases and executables</div>
    </div>    
    </a>
</li>
<li>
    <a class="prksLink" href="Blast.cgi?CMD=Web&amp;PAGE_TYPE=BlastDocs&amp;DOC_TYPE=DeveloperInfo" title="Develop Blast">
    <div id="dev-bg-img" class="bg-img">        
         <div class="prksTitle">Use BLAST API</div>
         <div class="prksDescr">Call BLAST from your application</div>
    </div>    
    </a>
</li>
<li>
    <a class="prksLink" href="Blast.cgi?CMD=Web&amp;PAGE_TYPE=BlastDocs&amp;DOC_TYPE=CloudBlast" title="Cloud Blast">
    <div id="cloud-bg-img" class="bg-img">        
         <div class="prksTitle">Use BLAST in the cloud</div>
         <div class="prksDescr">Start an instance at a cloud provider</div>
    </div>    
    </a>
</li>
</ul>
</div>
            <div class="links section">
<h2>Specialized searches</h2>
<ul>
<li class="spFirst"><a class="spLink" href="https://blast.ncbi.nlm.nih.gov/smartblast/?LINK_LOC=BlastHomeLink" title="SmartBlast">
    <div class="specblast-bg-img">        
         <div class="spTitle">SmartBLAST</div>
         <div class="spDesc">Find proteins highly similar to your query</div>
    </div>    
    </a>
</li>
<li><a class="spLink" <a href="https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=BlastHome" title="Find primers">
    <div class="specblast-bg-img">        
         <div class="spTitle">Primer-BLAST</div>
         <div class="spDesc">Design primers specific to your PCR template</div>
    </div>    
    </a>
</li>

<li><a class="spLink" href="Blast.cgi?PAGE_TYPE=BlastSearch&amp;PROG_DEF=blastn&amp;BLAST_PROG_DEF=blastn&amp;BLAST_SPEC=GlobalAln&amp;LINK_LOC=BlastHomeLink" title="Global alignment">
    <div class="specblast-bg-img">        
         <div class="spTitle">Global Align</div>
         <div class="spDesc">Compare two sequences across their entire span (Needleman-Wunsch)</div>
    </div>    
    </a>
</li>


<li><a class="spLink"href="https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi" title="Search conserved domains on a protein">
    <div class="specblast-bg-img">        
         <div class="spTitle">CD-search</div>
         <div class="spDesc">Find conserved domains in your sequence</div>
    </div>    
    </a>
</li>
<li class="spFirst"><a class="spLink"href="https://www.ncbi.nlm.nih.gov/igblast/" title="Search immunoglobulins">
    <div class="specblast-bg-img">        
         <div class="spTitle">IgBLAST</div>
         <div class="spDesc">Search immunoglobulins and T cell receptor sequences</div>
    </div>    
    </a>
</li>
<li>
<a class="spLink"href="https://www.ncbi.nlm.nih.gov/tools/vecscreen/" title="Screen a Sequence Using VecScreen">
    <div class="specblast-bg-img">        
         <div class="spTitle">VecScreen</div>
         <div class="spDesc">Search sequences for vector contamination</div>
    </div>    
    </a>
</li>
<li><a class="spLink"href="https://www.ncbi.nlm.nih.gov/Structure/lexington/lexington.cgi?cmd=rps" title="Find a protein sequence's domain architecture">
    <div class="specblast-bg-img">        
         <div class="spTitle">CDART</div>
         <div class="spDesc">Find sequences with similar conserved domain architecture</div>
    </div>    
    </a>
</li>
<li>
<a class="spLink"title="Compute a multiple protein sequence alignment" href="https://www.ncbi.nlm.nih.gov/tools/cobalt/cobalt.cgi?LINK_LOC=BlastHomeLink">
    <div class="specblast-bg-img">        
         <div class="spTitle">Multiple Alignment</div>
         <div class="spDesc">Align sequences using domain and protein constraints</div>
    </div>    
    </a>
</li>
<li class="spFirst">
<a class="spLink"title="Multiple Ordered Locus Estimates Tool" href="https://blast.ncbi.nlm.nih.gov/moleblast/moleblast.cgi">
    <div class="specblast-bg-img">        
         <div class="spTitle">MOLE-BLAST</div>
         <div class="spDesc">Establish taxonomy for uncultured or environmental sequences</div>
    </div>    
    </a>
</li>
</ul>
</div>                    
            
        </div><!--/content-->
     </div><!--/content-wrap-->        
</div><!--/wrap-->
 <footer>
      <section class="icon-section">
        <div id="icon-section-header" class="icon-section_header">Follow NCBI</div>
        <div class="grid-container container">
          <div class="icon-section_container">
            <a class="footer-icon" id="footer_twitter" href="https://twitter.com/ncbi" aria-label="Twitter"><svg data-name="Layer 1" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 300 300">
                <defs>
                  <style>
                    .cls-11 {
                      fill: #737373;
                    }
                  </style>
                </defs>
                <title>Twitter</title>
                <path class="cls-11" d="M250.11,105.48c-7,3.14-13,3.25-19.27.14,8.12-4.86,8.49-8.27,11.43-17.46a78.8,78.8,0,0,1-25,9.55,39.35,39.35,0,0,0-67,35.85,111.6,111.6,0,0,1-81-41.08A39.37,39.37,0,0,0,81.47,145a39.08,39.08,0,0,1-17.8-4.92c0,.17,0,.33,0,.5a39.32,39.32,0,0,0,31.53,38.54,39.26,39.26,0,0,1-17.75.68,39.37,39.37,0,0,0,36.72,27.3A79.07,79.07,0,0,1,56,223.34,111.31,111.31,0,0,0,116.22,241c72.3,0,111.83-59.9,111.83-111.84,0-1.71,0-3.4-.1-5.09C235.62,118.54,244.84,113.37,250.11,105.48Z">
                </path>
              </svg></a>
            <a class="footer-icon" id="footer_facebook" href="https://www.facebook.com/ncbi.nlm" aria-label="Facebook"><svg data-name="Layer 1" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 300 300">
                <title>Facebook</title>
                <path class="cls-11" d="M210.5,115.12H171.74V97.82c0-8.14,5.39-10,9.19-10h27.14V52l-39.32-.12c-35.66,0-42.42,26.68-42.42,43.77v19.48H99.09v36.32h27.24v109h45.41v-109h35Z">
                </path>
              </svg></a>
            <a class="footer-icon" id="footer_linkedin" href="https://www.linkedin.com/company/ncbinlm" aria-label="LinkedIn"><svg data-name="Layer 1" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 300 300">
                <title>LinkedIn</title>
                <path class="cls-11" d="M101.64,243.37H57.79v-114h43.85Zm-22-131.54h-.26c-13.25,0-21.82-10.36-21.82-21.76,0-11.65,8.84-21.15,22.33-21.15S101.7,78.72,102,90.38C102,101.77,93.4,111.83,79.63,111.83Zm100.93,52.61A17.54,17.54,0,0,0,163,182v61.39H119.18s.51-105.23,0-114H163v13a54.33,54.33,0,0,1,34.54-12.66c26,0,44.39,18.8,44.39,55.29v58.35H198.1V182A17.54,17.54,0,0,0,180.56,164.44Z">
                </path>
              </svg></a>
            <a class="footer-icon" id="footer_github" href="https://github.com/ncbi" aria-label="GitHub"><svg data-name="Layer 1" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 300 300">
                <defs>
                  <style>
                    .cls-11,
                    .cls-12 {
                      fill: #737373;
                    }

                    .cls-11 {
                      fill-rule: evenodd;
                    }
                  </style>
                </defs>
                <title>GitHub</title>
                <path class="cls-11" d="M151.36,47.28a105.76,105.76,0,0,0-33.43,206.1c5.28,1,7.22-2.3,7.22-5.09,0-2.52-.09-10.85-.14-19.69-29.42,6.4-35.63-12.48-35.63-12.48-4.81-12.22-11.74-15.47-11.74-15.47-9.59-6.56.73-6.43.73-6.43,10.61.75,16.21,10.9,16.21,10.9,9.43,16.17,24.73,11.49,30.77,8.79,1-6.83,3.69-11.5,6.71-14.14C108.57,197.1,83.88,188,83.88,147.51a40.92,40.92,0,0,1,10.9-28.39c-1.1-2.66-4.72-13.42,1-28,0,0,8.88-2.84,29.09,10.84a100.26,100.26,0,0,1,53,0C198,88.3,206.9,91.14,206.9,91.14c5.76,14.56,2.14,25.32,1,28a40.87,40.87,0,0,1,10.89,28.39c0,40.62-24.74,49.56-48.29,52.18,3.79,3.28,7.17,9.71,7.17,19.58,0,14.15-.12,25.54-.12,29,0,2.82,1.9,6.11,7.26,5.07A105.76,105.76,0,0,0,151.36,47.28Z">
                </path>
                <path class="cls-12" d="M85.66,199.12c-.23.52-1.06.68-1.81.32s-1.2-1.06-.95-1.59,1.06-.69,1.82-.33,1.21,1.07.94,1.6Zm-1.3-1">
                </path>
                <path class="cls-12" d="M90,203.89c-.51.47-1.49.25-2.16-.49a1.61,1.61,0,0,1-.31-2.19c.52-.47,1.47-.25,2.17.49s.82,1.72.3,2.19Zm-1-1.08">
                </path>
                <path class="cls-12" d="M94.12,210c-.65.46-1.71,0-2.37-.91s-.64-2.07,0-2.52,1.7,0,2.36.89.65,2.08,0,2.54Zm0,0"></path>
                <path class="cls-12" d="M99.83,215.87c-.58.64-1.82.47-2.72-.41s-1.18-2.06-.6-2.7,1.83-.46,2.74.41,1.2,2.07.58,2.7Zm0,0">
                </path>
                <path class="cls-12" d="M107.71,219.29c-.26.82-1.45,1.2-2.64.85s-2-1.34-1.74-2.17,1.44-1.23,2.65-.85,2,1.32,1.73,2.17Zm0,0">
                </path>
                <path class="cls-12" d="M116.36,219.92c0,.87-1,1.59-2.24,1.61s-2.29-.68-2.3-1.54,1-1.59,2.26-1.61,2.28.67,2.28,1.54Zm0,0">
                </path>
                <path class="cls-12" d="M124.42,218.55c.15.85-.73,1.72-2,1.95s-2.37-.3-2.52-1.14.73-1.75,2-2,2.37.29,2.53,1.16Zm0,0"></path>
              </svg></a>
            <a class="footer-icon" id="footer_blog" href="https://ncbiinsights.ncbi.nlm.nih.gov/" aria-label="Blog">
              <svg id="Layer_1" data-name="Layer 1" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 40 40"><defs><style>.cls-1{fill:#737373;}</style></defs><path class="cls-1" d="M14,30a4,4,0,1,1-4-4,4,4,0,0,1,4,4Zm11,3A19,19,0,0,0,7.05,15a1,1,0,0,0-1,1v3a1,1,0,0,0,.93,1A14,14,0,0,1,20,33.07,1,1,0,0,0,21,34h3a1,1,0,0,0,1-1Zm9,0A28,28,0,0,0,7,6,1,1,0,0,0,6,7v3a1,1,0,0,0,1,1A23,23,0,0,1,29,33a1,1,0,0,0,1,1h3A1,1,0,0,0,34,33Z"></path></svg>
            </a>
          </div>
        </div>
      </section>

      <section class="container-fluid bg-primary">
        <div class="container pt-5">
          <div class="row mt-3">
            <div class="col-lg-3 col-12">
              <p><a class="text-white" href="https://www.nlm.nih.gov/socialmedia/index.html">Connect with NLM</a></p>
              <ul class="list-inline social_media">
                <li class="list-inline-item"><a href="https://twitter.com/NLM_NIH" aria-label="Twitter" target="_blank" rel="noopener noreferrer"><svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" viewBox="0 0 249 249" style="enable-background:new 0 0 249 249;" xml:space="preserve">
                      <style type="text/css">
                        .st20 {
                          fill: #FFFFFF;
                        }

                        .st30 {
                          fill: none;
                          stroke: #FFFFFF;
                          stroke-width: 8;
                          stroke-miterlimit: 10;
                        }
                      </style>
                      <title>SM-Twitter</title>
                      <g>
                        <g>
                          <g>
                            <path class="st20" d="M192.9,88.1c-5,2.2-9.2,2.3-13.6,0.1c5.7-3.4,6-5.8,8.1-12.3c-5.4,3.2-11.4,5.5-17.6,6.7
                                                c-10.5-11.2-28.1-11.7-39.2-1.2c-7.2,6.8-10.2,16.9-8,26.5c-22.3-1.1-43.1-11.7-57.2-29C58,91.6,61.8,107.9,74,116
                                                c-4.4-0.1-8.7-1.3-12.6-3.4c0,0.1,0,0.2,0,0.4c0,13.2,9.3,24.6,22.3,27.2c-4.1,1.1-8.4,1.3-12.5,0.5c3.6,11.3,14,19,25.9,19.3
                                                c-11.6,9.1-26.4,13.2-41.1,11.5c12.7,8.1,27.4,12.5,42.5,12.5c51,0,78.9-42.2,78.9-78.9c0-1.2,0-2.4-0.1-3.6
                                                C182.7,97.4,189.2,93.7,192.9,88.1z"></path>
                          </g>
                        </g>
                        <circle class="st30" cx="124.4" cy="128.8" r="108.2"></circle>
                      </g>
                    </svg></a></li>
                <li class="list-inline-item"><a href="https://www.facebook.com/nationallibraryofmedicine" aria-label="Facebook" rel="noopener noreferrer" target="_blank">
                    <svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" viewBox="0 0 249 249" style="enable-background:new 0 0 249 249;" xml:space="preserve">
                      <style type="text/css">
                        .st10 {
                          fill: #FFFFFF;
                        }

                        .st110 {
                          fill: none;
                          stroke: #FFFFFF;
                          stroke-width: 8;
                          stroke-miterlimit: 10;
                        }
                      </style>
                      <title>SM-Facebook</title>
                      <g>
                        <g>
                          <path class="st10" d="M159,99.1h-24V88.4c0-5,3.3-6.2,5.7-6.2h16.8V60l-24.4-0.1c-22.1,0-26.2,16.5-26.2,27.1v12.1H90v22.5h16.9
                                                      v67.5H135v-67.5h21.7L159,99.1z"></path>
                        </g>
                      </g>
                      <circle class="st110" cx="123.6" cy="123.2" r="108.2"></circle>
                    </svg>
                  </a></li>
                <li class="list-inline-item"><a href="https://www.youtube.com/user/NLMNIH" aria-label="Youtube" target="_blank" rel="noopener noreferrer"><svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" viewBox="0 0 249 249" style="enable-background:new 0 0 249 249;" xml:space="preserve">
                      <title>SM-Youtube</title>
                      <style type="text/css">
                        .st4 {
                          fill: none;
                          stroke: #FFFFFF;
                          stroke-width: 8;
                          stroke-miterlimit: 10;
                        }

                        .st5 {
                          fill: #FFFFFF;
                        }
                      </style>
                      <circle class="st4" cx="124.2" cy="123.4" r="108.2"></circle>
                      <g transform="translate(0,-952.36218)">
                        <path class="st5" d="M88.4,1037.4c-10.4,0-18.7,8.3-18.7,18.7v40.1c0,10.4,8.3,18.7,18.7,18.7h72.1c10.4,0,18.7-8.3,18.7-18.7
                                            v-40.1c0-10.4-8.3-18.7-18.7-18.7H88.4z M115.2,1058.8l29.4,17.4l-29.4,17.4V1058.8z"></path>
                      </g>
                    </svg></a></li>
              </ul>
            </div>
            <div class="col-lg-3 col-12">
              <p class="address_footer text-white">National Library of Medicine<br>
                <a href="https://www.google.com/maps/place/8600+Rockville+Pike,+Bethesda,+MD+20894/@38.9959508,-77.101021,17z/data=!3m1!4b1!4m5!3m4!1s0x89b7c95e25765ddb:0x19156f88b27635b8!8m2!3d38.9959508!4d-77.0988323" class="text-white" target="_blank" rel="noopener noreferrer">8600 Rockville Pike<br>
                  Bethesda, MD 20894</a></p>
            </div>
            <div class="col-lg-3 col-12 centered-lg">
              <p><a href="https://www.nlm.nih.gov/web_policies.html" class="text-white">Web Policies</a><br>
                <a href="https://www.nih.gov/institutes-nih/nih-office-director/office-communications-public-liaison/freedom-information-act-office" class="text-white">FOIA</a><br>
                <a href="https://www.hhs.gov/vulnerability-disclosure-policy/index.html" class="text-white" id="vdp">HHS Vulnerability Disclosure</a></p>
            </div>
            <div class="col-lg-3 col-12 centered-lg">
              <p><a class="supportLink text-white" href="https://support.nlm.nih.gov/">Help</a><br>
                <a href="https://www.nlm.nih.gov/accessibility.html" class="text-white">Accessibility</a><br>
                <a href="https://www.nlm.nih.gov/careers/careers.html" class="text-white">Careers</a></p>
            </div>
          </div>
          <div class="row">
            <div class="col-lg-12 centered-lg">
              <nav class="bottom-links">
                <ul class="mt-3">
                  <li>
                    <a class="text-white" href="//www.nlm.nih.gov/">NLM</a>
                  </li>
                  <li>
                    <a class="text-white" href="https://www.nih.gov/">NIH</a>
                  </li>
                  <li>
                    <a class="text-white" href="https://www.hhs.gov/">HHS</a>
                  </li>
                  <li>
                    <a class="text-white" href="https://www.usa.gov/">USA.gov</a>
                  </li>
                </ul>
              </nav>
            </div>
          </div>
        </div>
      </section>
    </footer>
<script type="text/javascript" src="js/nwds.js"></script>
<script type="text/javascript" src="js/ncbipopup.js"></script>
<script type="text/javascript" src="js/headerNew.js"></script>
<script src="js/uswds.min.js"></script>
<script type="text/javascript" src="/portal/portal3rc.fcgi/rlib/js/InstrumentOmnitureBaseJS/InstrumentNCBIBaseJS/InstrumentPageStarterJS.js"></script>
<!--
<script type="text/javascript" src="://www.ncbi.nlm.nih.gov/portal/portal3rc.fcgi/supportedbrowsers/js/nonportalbc_min.js"></script>
<script type="text/javascript">$("#browsers_ajax").browserCheck();</script>
-->
</body>
</html>
