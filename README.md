# [Unsupervised Distance Learning Framework](http://www.ic.unicamp.br/~dcarlos/UDLF/index.html)

[Access the official software webpage >>](http://www.ic.unicamp.br/~dcarlos/UDLF/index.html)

**Authors:** [Lucas Pascotti Valem](http://www.lucasvalem.com) and [Daniel Carlos Guimarães Pedronette](http://www.ic.unicamp.br/~dcarlos/)

Dept. of Statistic, Applied Math. and Computing, Universidade Estadual Paulista ([UNESP](http://www.rc.unesp.br/)), Rio Claro, Brazil

----------------------
* [Overview](#overview)
* [Get Started](#get-started)
* [Binaries](#binaries)
* [Compilation](#compilation)
* [Execution](#execution)
* [Documentation](#documentation)
* [Contributing](#contributing)
* [Cite](#cite)
* [Contact](#contact)
* [Acknowledgments](#acknowledgments)
* [License](#license)

## Overview
A framework of unsupervised distance learning methods for image and multimedia retrieval tasks.
Currently, eleven different [unsupervised learning methods](https://github.com/UDLF/UDLF/wiki/Methods) are implemented
([RFE](https://doi.org/10.1109/TIP.2023.3268868),
[RDPAC](https://doi.org/10.3390/jimaging7030049),
[BFSTree](https://doi.org/10.1016/j.patcog.2020.107666),
[LHRR](http://doi.org/10.1109/TIP.2019.2920526),
[ContextRR](http://dl.acm.org/citation.cfm?id=1948207.1948291),
[Correlation Graph](http://dx.doi.org/10.1016/j.neucom.2016.03.081),
[CPRR](http://dx.doi.org/10.1109/SIBGRAPI.2016.042),
[Rk Graph Dist.](http://dx.doi.org/10.1016/j.patrec.2016.05.021),
[ReckNNGraph](http://dx.doi.org/10.1016/j.imavis.2013.12.009),
[RL-Recom](http://dx.doi.org/10.1145/2671188.2749336),
and [RL-Sim*](http://dx.doi.org/10.1145/2671188.2749335)).

## Get Started
An easy guide for your first use can be found in the [software official webpage](http://www.ic.unicamp.br/%7Edcarlos/UDLF/getStarted.html).

## Binaries
Binaries are available for download in the [release page](https://github.com/UDLF/UDLF/releases).

## Compilation
This project can be compiled by any C++ compiler that supports the C++2014 standard. There is a Makefile that can be used to compile the code. A executable called `udlf` will be generate inside the **bin/** directory.

## Execution
The executable is called in the terminal:

- **Linux and MacOS:** `./udlf [config.ini]`

- **Windows:**  `call udlf.exe [config.ini]`

**NOTE:** The binary must be executed inside the **bin/** directory.

The [configuration file](https://github.com/UDLF/UDLF/wiki/Configuration) specifies everything about the execution:
the desired task, method being used, dataset information, [input files](https://github.com/UDLF/UDLF/wiki/File-Formats),
[output files](https://github.com/UDLF/UDLF/wiki/File-Formats),
[evaluation settings](https://github.com/UDLF/UDLF/wiki/Evaluation),
and other details.
When the binary is executed, it searchs for a `config.ini` file in its current directory. The user can also specify a different
configuration file that can be passed as a parameter: `./udlf my_conf.ini.` The software considers only a single configuration file per execution.

**NOTE:** Complete examples of input files for distinct datasets are available [here](https://github.com/UDLF/Datasets).

After the execution, a `log.txt` is generated:

```
 - GENERAL INFORMATION -
 --------------------------------------
 Task:             UDL
 Method:           CPRR
 Dataset Size:     1400
 Image List File:  desc/lists/mpeg7.txt
 Image Class File: desc/classes/mpeg7.txt
 Input File:       desc/matrices/mpeg7/cfd.txt
 Input Format:     MATRIX DIST
 Output File:      output/output
 Output Format:    RK ALL
 --------------------------------------
 - METHOD PARAMETERS -
 --------------------------------------
 PARAM_CPRR_K = 20
 PARAM_CPRR_L = 400
 PARAM_CPRR_T = 2
 --------------------------------------
 - EVALUATION RESULTS -
 --------------------------------------
 * Efficiency: Total Time of the Algorithm Execution: 0.0438 s
 * Effectiveness:
 Before:
	 P@20		0.7559
	 Recall@40	0.8444
	 MAP		0.8064
 After:
	 P@20		0.8979
	 Recall@40	0.9477
	 MAP		0.9215
 Relative Gains:
	 P@20		18.7866%
	 Recall@40	12.2404%
	 MAP		14.2707%
 --------------------------------------
 Log generated at 2017/1/26 16:37:24
```

The results can be exported in [different formats](https://github.com/UDLF/UDLF/wiki/File-Formats).
Below you can see some examples of ranked lists that were exported as a *html* page.
The query images are presented in green borders and wrong results in red borders.
The first line represents the original retrieval results and the second line, the results after the algorithm execution.

![corel5k](https://github.com/UDLF/UDLF/blob/master/visual_examples/corel5k.png)

![mpeg7](https://github.com/UDLF/UDLF/blob/master/visual_examples/mpeg7.png)

![oxford17flowers](https://github.com/UDLF/UDLF/blob/master/visual_examples/oxford17flowers.png)

![soccer](https://github.com/UDLF/UDLF/blob/master/visual_examples/soccer.png)

**NOTE:** The above examples consider the datasets
[Corel5k](http://www.ci.gxnu.edu.cn/cbir/Dataset.aspx),
[MPEG-7](http://www.dabi.temple.edu/~shape/MPEG7/dataset.html),
[Oxford17Flowers](http://www.robots.ox.ac.uk/~vgg/data/flowers/), and
[Soccer](http://lear.inrialpes.fr/people/vandeweijer/data.html);
respectively.

## Documentation
The documentation is available in the [software wiki.](https://github.com/UDLF/UDLF/wiki)

## Contributing
We appreciate suggestions, ideas and contributions.
If you want to contribute, feel free to [contact us.](#contact)
Github pull requests should be avoided because they are not part of our review process.
To report small bugs, you can use the [issue tracker](https://github.com/UDLF/UDLF/issues) provided by GitHub.

## Cite
If you use this software, please cite

 ```latex
@inproceedings{Valem:2017:UDL:3078971.3079017,
	author = {Valem, Lucas Pascotti and Pedronette, Daniel Carlos Guimar\~{a}es},
	title = {An Unsupervised Distance Learning Framework for Multimedia Retrieval},
	booktitle = {Proceedings of the 2017 ACM on International Conference on Multimedia Retrieval},
	series = {ICMR '17},
	year = {2017},
	isbn = {978-1-4503-4701-3},
	location = {Bucharest, Romania},
	pages = {107--111},
	numpages = {5},
	url = {http://doi.acm.org/10.1145/3078971.3079017},
	doi = {10.1145/3078971.3079017},
	acmid = {3079017},
	publisher = {ACM},
	address = {New York, NY, USA},
}
```

## Contact
**Lucas Pascotti Valem**: `lucaspascottivalem@gmail.com` or `lucas.valem@unesp.br`

**Daniel Carlos Guimarães Pedronette**: `daniel.pedronette@unesp.br`

## Acknowledgments
The authors are grateful to São Paulo Research Foundation - [FAPESP](http://www.fapesp.br/en/) (grants 2013/08645-0, and 2014/04220-8).

## License
This project is licensed under GPLv2. See [details.](https://github.com/UDLF/UDLF/blob/master/LICENSE)
