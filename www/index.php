<!doctype html>
<html>
<head>
<meta charset="UTF-8">

<link rel="apple-touch-icon" sizes="180x180" href="/apple-touch-icon.png">
<link rel="icon" type="image/png" href="/favicon-32x32.png" sizes="32x32">
<link rel="icon" type="image/png" href="/favicon-16x16.png" sizes="16x16">
<link rel="manifest" href="/manifest.json">
<link rel="mask-icon" href="/safari-pinned-tab.svg" color="#5bbad5">
<meta name="theme-color" content="#ffffff">

<link rel="stylesheet" href="style480.css" media="screen and (min-device-width: 320px) and (max-device-width: 667px)" />
<link rel="stylesheet" href="style780.css" media="screen and (min-device-width: 668px) and (max-device-width: 959px)" />
<link rel="stylesheet" href="style960.css" media="screen and (min-device-width: 960px)" />

<title>Project <i>lefko3</i></title>
<style>

</style>
</head>

<body>

<div id = "header"><br>

<picture> 
    <source media = "(max-device-width: 480pt)" srcset = "images/header3a.png" class = "centerImage">
	<source media = "(max-device-width: 780pt)" srcset = "images/header3c.png" class = "centerImage">
    <img src = "images/header3d.png" class = "centerImage" alt="banner image">
</picture>

<nav>
    <a href = "http://www.sheffersonlab.com/index.html"><picture>  
        <source media = "(max-device-width: 667px)" srcset = "images/home3.png" class = "NavImg">
        <img src = "images/home1.png" alt="Home" class = "NavImg">
    </picture></a>
	<a href = "http://www.sheffersonlab.com/about.html"><picture>  
        <source media = "(max-device-width: 667px)" srcset = "images/about3.png" class = "NavImg">
        <img src = "images/about1.png" alt="About" class = "NavImg">
    </picture></a>
	<a href = "http://www.sheffersonlab.com/news.html"><picture>
        <source media = "(max-device-width: 667px)" srcset = "images/news3.png" class = "NavImg">
        <img src = "images/news1.png" alt="News" class = "NavImg">
	</picture></a>
 	<a href = "http://www.sheffersonlab.com/researchp.html"><picture>
        <source media = "(max-device-width: 667px)" srcset = "images/research3.png" class = "NavImg">
        <img src = "images/research1.png" alt="Research" class = "NavImg">
    </picture></a>
	<a href = "http://www.sheffersonlab.com/labpubs.html"><picture>
        <source media = "(max-device-width: 667px)" srcset = "images/pubs3.png" class = "NavImg">
        <img src = "images/pubs1.png" alt="Publications" class = "NavImg">
	</picture></a>
	<a href = "http://www.sheffersonlab.com/people.html"><picture>
        <source media = "(max-device-width: 667px)" srcset = "images/peops3.png" class = "NavImg">
        <img src = "images/peops1.png" alt="People" class = "NavImg">
	</picture></a>
	<a href = "http://www.sheffersonlab.com/projects.html"><picture>
        <source media = "(max-device-width: 667px)" srcset = "images/events4.png" class = "NavImg">
        <img src = "images/events2.png" alt="Events" class = "NavImg">
    </picture></a>
	<a href = "http://www.sheffersonlab.com/classes.html"><picture>
        <source media = "(max-device-width: 667px)" srcset = "images/classes3.png" class = "NavImg">
        <img src = "images/classes1.png" alt="Classes" class = "NavImg">
	</picture></a>
</nav>
    
<picture> 
    <source media = "(max-device-width: 480px)" srcset = "images/grousepointthin.png" class = "centerImage">
	<source media = "(max-device-width: 780px)" srcset = "images/IMG_1571.p.780.jpg" class = "centerImage">
    <img src = "images/banner3.png" class = "centerImage" alt="banner image" width = "950px">
</picture>


</div>

<div id = "main">
	<ul>
		<h1>Summary of R Package  <i>lefko3</i></h1>
		<p><i>Lefko3</i> is an R package designed to provide a complete working environment to create and analyze size-classified matrix projection models (MPMs), including Lefkovitch MPMs, integral projection models (IPMs), and age-by-stage MPMs, in both raw and function-based formats. It handles the entire workflow in MPM development, from formatting the demographic dataset, to building the life cycle graph, to parameterizing vital rate models, to contructing the MPM, and finally to analysis. It was originally designed to produce historical and ahsitorical matrix projection models (MPMs), where historical MPMs are large matrices that incorporate 3 consecutive time steps in each matrix. These matrices allow the history of the individual to be more completely modeled within population dynamics models, which is important given that MPMs typically assume no influence of differences in individual history on population dynamics. However, it now handles all forms of MPM.</p>
		<p>Package <i>lefko3</i> works on R 3.6.3 or greater on all systems that support R, including Windows, MacOS, Linux, and Unix.</p>
		<p>The project is hosted on R-Forge, with the most recent stable version of <i>lefko3</i> available for download through CRAN, with documentation available <a href = "https://cran.r-project.org/package=lefko3">here</a>. The most recent development version is available through R-Forge <a href = "https://r-forge.r-project.org/R/?group_id=2427">here</a>. Please also see our paper and cite our paper describing <i>lefko3</i> in <a href = "https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13526">Methods in Ecology and Evolution</a>.</p>
		<br>
		<p><b>Current version:</b> 3.6.0 (available from July 2021).</p>
		<p><b>New in this version:</b> Sped up analysis functions through better memory management.</p>
		<p><b>Next version:</b> 3.7.0 (expected late-August 2021). </p>
		<p><b>Coming in next version:</b>  LTRE and sLTRE functions. New functions to create and edit lefkoMat objects.</p>
		<p><b>On the drawing board:</b> Bootstrapping functions for vertical datasets. Parallelization functions.</p>
		<p><b>Tutorial materials:</b>
			<ol>
				<li>Tutorial on <a href = "https://cran.r-project.org/web/packages/lefko3/vignettes/Chapter1.html">theory, methods, and concepts</a> (included in package on CRAN)</li>
				<li>Tutorial example of <a href = "https://cran.r-project.org/web/packages/lefko3/vignettes/Chapter2.html">raw historical and ahistorical MPM analysis using <i>Lathyrus vernus</i></a> (included in package on CRAN)</li>
				<li>Tutorial example of <a href = "https://cran.r-project.org/web/packages/lefko3/vignettes/Chapter5.html">function-based historical and ahistorical MPM analysis using <i>Cypripedium candidum</i></a> (included in package on CRAN)</li>
				<li>Tutorial example of <a href = "https://cran.r-project.org/web/packages/lefko3/vignettes/Chapter6.html">historical and ahistorical IPM analysis using <i>Lathyrus vernus</i></a> (included in package on CRAN)</li>
				<li>Tutorial example of <a href = "https://cran.r-project.org/web/packages/lefko3/vignettes/Chapter8.html">ahistorical age x stage MPM analysis using <i>Lathyrus vernus</i></a> (included in package on CRAN)</li>
				<li>Video tutorial: <a href = "https://youtu.be/0MtIcXgarEw">Intro to <i>lefko3</i> with a short analysis of <i>Cypripedium candidum</i></a></li>
			</ol>
		</p>
	<p><b>Current workshop schedule:</b></p>
		<ol>
			<li>Power Demography Workshop - an Intro to <i>lefko3</i>. 2 Sept 2021 at 4am GMT (1pm Tokyo time).</li>
			<li>Expected workshop in mid September for North / South American users.</li>
			<li>Expected workshop in mid September for European / African / Central Asian users.</li>
		</ol>
	
	<p>Please contact the lab for further details.</p>
		<br>
		<br>
		<br>
	</ul>

	<br>

</div>
<br>

<div id = "rightcolumn">
	<a class="twitter-timeline" data-width="90%" data-height="400pt" href="https://twitter.com/Richshefferson" data-chrome="noheader nofooter noborders noscrollbar">Tweets by Richshefferson</a> <script async src="//platform.twitter.com/widgets.js" charset="utf-8" float = "right"></script>
	<a href="https://twitter.com/intent/tweet?screen_name=Richshefferson" class="twitter-mention-button" data-show-count="false">Tweet to @Richshefferson</a><script async src="//platform.twitter.com/widgets.js" charset="utf-8"></script>
	<br>
</div>

<div id = "footer">
	<hr>
   	<br>
	<p>The Shefferson lab is recruiting students and post-docs!</p>
	<p>Interested parties may contact us through the <a href = "http://gpes.c.u-tokyo.ac.jp/">GPES</a> and <a href = "http://system.c.u-tokyo.ac.jp/p-graduate/guide.html">広域システム</a> applicant pages.</p>
	<br>
	<br>
	<p>All content copyright &copy; 2021, Richard P. Shefferson<br>
	Lady's slipper painting courtesy of Carmen Kraus</p>
</div>

</body>
</html>
