#' KitchenhamEtAl.CorrelationsAmongParticipants.Madeyski10 data
#'
#' Data illustrate correlations between results from individual participants in cross-over
#' experiment P2007 (Smell and Library) conducted by Madeyski, see:
#' [1] Lech Madeyski, Test-Driven Development: An Empirical Evaluation of Agile Practice.
#' (Heidelberg, London, New York): Springer, 2010. Foreword by Prof. Claes Wohlin.
#
#' If you use this data set please cite:
#' [1] Lech Madeyski, Test-Driven Development: An Empirical Evaluation of Agile Practice.
#' (Heidelberg, London, New York): Springer, 2010. Foreword by Prof. Claes Wohlin.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format `KitchenhamEtAl.CorrelationsAmongParticipants.Madeyski10`: a data frame with 45 rows
#' and 10 variables:
#' \describe{
#' \item{ExperimentID}{<fct>| ExperimentID: This experiment is the only cross-over experiment in
#' the family of TDD and Pair-Programming experiments conducted by Madeyski, so all values in this
#' column are set to 'P2007'.}
#' \item{ParticipantID}{<fct> | Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct> | Experimental Sequence Group: A (TLSP-TFSP), B (TFSP-TLSP)}
#' \item{System}{<fct> | Software system to develop: Smell (a tool for identifying bad code smells
#' in Java source code through the use of a set of software metrics) or Library (a library
#' application)}
#' \item{Period}{<fct> | Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct> | Experimental Treatment: Test-First Solo Programming (TFSP) vs Test-Last
#' Solo Programming (TLSP)}
#' \item{PATP}{<dbl> | Dependent variable: Percentage of Acceptance Tests Passed}
#' \item{NATPPH}{<dbl> | Dependent variable: Number of Acceptance Tests Passed Per Hour}
#' \item{CBO}{<dbl> | Dependent variable: Mean value of Coupling Between Objects (CBO), see CK set
#' of metrics}
#' \item{WMC}{<dbl> | Dependent variable: Mean value of Weighted Number of Methods in Class (WMC),
#' see CK set of metrics}
#' \item{RFC}{<dbl> | Dependent variable: Mean value of Response For a Class (RFC), see CK set of
#' metrics}
#' \item{CrossOverID}{<fct> | Cross-Over Code. This experiment is a simple two-group cross-over
#' experiment with one cross-over code, so all values in this column are set to 'CO1'.
#' However, four-group experiments require a code to identify the linked sequence groups
#' (although that can be deduced from the system used in the first time period).
#' A crossover code is also essential for non-parametric analysis.}
#' }
#' @source \url{https://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Madeyski10
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Madeyski10"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello15EMSE data
#'
#' Data illustrate correlations between results from individual participants in cross-over
#' experiment usb2 conducted by Scanniello et al:
#' [1] G. Scanniello, A. Marcus, and D. Pascale, 'Link analysis algorithms for static concept
#' location: an empirical assessment', Empirical Software Engineering, vol. 20, no. 6,
#' pp. 1666–1720, 2015.
#' The goal of the experiment is to assess whether a new technique (implemented as an Eclipse
#' plug-in) for static concept location (proposed by the authors) supports users in identifying the
#' places in the code where changes are to be made.
#'
#' If you use this data set please cite:
#' [1] G. Scanniello, A. Marcus, and D. Pascale, 'Link analysis algorithms for static concept
#' location: an empirical assessment', Empirical Software Engineering, vol. 20, no. 6,
#' pp. 1666–1720, 2015.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format `KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello15EMSE`: a data frame with 48
#' rows and 10 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each experiment in the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A (CL-NOCL,Jedit-Atunes),
#' B (NOCL-CL,Atunes-Jedit), C(NOCL-CL,Jedit-Atunes), D(CL-NOCL,Atunes-Jedit)}
#' \item{System}{<fct>|Software systems used in the experiment: Jedit and Atunes}
#' \item{Treatment}{<fct>|Experimental Treatment: Use of Concept Location plug-in (CL) vs
#' no Concept Location plug-in (NOCL)}
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Correctness}{<int>|Dependent variable: 0, 1, 2, 3, 4. The participants are asked to
#' indicate a single change method for each of 4 bug reports. A change method is correctly
#' identified if that method is in the change set of the bug report.}
#' \item{Time}{<dbl>|Dependent variable: The total time [min] to accomplish concept location tasks,
#' i.e.,to identify (four) bugs given their reports}
#' \item{Efficiency}{<dbl>|Dependent variable: The participants’ efficiency in the execution of
#' concept location tasks. It is computed dividing correctness by time. }
#' \item{CrossOverID}{<fct>|Crossover category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. }
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello15EMSE
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello15EMSE"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM data
#'
#' Data illustrate correlations between results from individual participants in a family of four
#' cross-over experiments conducted by Scanniello et al:
#' [1] G.  Scanniello,  C. Gravino,  M. Genero, J.A. Cruz-Lemus, and  G. Tortora,  'On the Impact
#' of UML Analysis Models on Source-Code Comprehensibility and Modifiability', ACM Transactions on
#' Software Engineering and Methodlogy, vol. 23, no. 2, pp. 13:1-13:26, 2014
#' The family of experiments investigated whether the availability of analysis models in addition
#' to the source code made the code easier to understand and modify.
#' If you use this data set please cite:
#' [1] G.  G.  Scanniello,  C. Gravino,  M. Genero, J.A. Cruz-Lemus, and  G. Tortora, 'On the
#' Impact of UML Analysis Models on Source-Code Comprehensibility and Modifiability', ACM
#' Transactions on Software Engineering and Methodology, vol. 23, no. 2, pp. 13:1-13:26, 2014
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#' @format `KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM`: a data frame with 172
#' rows and 9 variables:
#' \describe{
#' \item{ExperimentID}{<fct> | ExperimentID: A unique identifier for each experiment in the data
#' set.}
#' \item{ParticipantID}{<fct> | Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{Treatment}{<fct> | Experimental Treatment: AM an Analysis Model with source code (AM) vs
#' Source Code only (SC)}
#' \item{SequenceGroup}{<fct> | Experimental Sequence Group: A (AM-SC,S1-S2), B (SC-AM,S1-S2),
#' C(AM-SC,S2-S1), D(SC-AM,S2-S1)}
#' \item{System}{<fct> | Software systems used in the experiment: S1 A system to sell and manage
#' CDs/DVDs in a music shop, S2 A system to book and by theater tickets.}
#' \item{Comprehension}{<dbl> | Dependent variable: The comprehension level the software engineer
#' achieved based on the F-measure }
#' \item{Modification}{<dbl> | Dependent variable: The modifiability level the software engineer
#' achieved based on the F-measure}
#' \item{Period}{<fct> | Time period of the cross-over experiment: 1 or 2}
#' \item{CrossOverID}{<fct> | CrossOver category: For 4 group the crossover category specifies the
#' matching pairs of sequence groups, CO1 and CO2. For 2 group crossover, the category is set to
#' CO1 only}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Torchiano17JVLC data
#'
#' Data illustrate correlations between results from individual participants in a family of three
#' cross-over experiments conducted by Torchiano et al:
#' [1] M. Torchiano, G. Scanniello, F. Ricca, G. Reggio, and M. Leotta,
#' 'Do UML object diagrams affect design comprehensibility? Results
#' from a family of four controlled experiments.' Journal of Visual Languages
#' and  Computing, vol. 41, pp. 10–21, 2017.
#' Although the paper reports four experiment, we only have data from three of those experiments.
#' The experiments assess whether the comprehensibility of UML specifications improve when the
#' software documents include UML object diagrams as well as the standard UML class diagrams.
#' If you use this data set please cite:
#' [1] M. Torchiano, G. Scanniello, F. Ricca, G. Reggio, and M. Leotta,
#' 'Do UML object diagrams affect design comprehensibility? Results
#' from a family of four controlled experiments.' Journal of Visual Languages
#' and  Computing, vol. 41, pp. 10–21, 2017.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format A data frame with 214 rows and 8 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each of the three experiments
#' in the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A , B , C, D}
#' \item{System}{<fct>|Software systems used in the experiment: File System manager (FS) for
#' folders, files, links. Roads system (R) handles maps made up of cities connected by means of
#' roads. Train (T) a system to manage timetables, trains, and paths. Catalogue system (C).
#' It collects category of items (e.g., cars) and items (e.g., car models) based on a set of
#' features (e.g., number of  doors). In PoliTo2, only FS and T were administered to the
#' participants, while in UniBas1 and UniGe1 all the four experimental objects were used. }
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment: Object Diagram (OD) vs No Object Diagram (NoOD)}
#' \item{Comprehension}{<dbl>|Dependent variable: The comprehension level the software engineers.
#' For PoliTo2 Comprehension was based on answering a set of 4 questions, for UniBas and UniGe
#' comprehension was measured using the F metric. }
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. For 2 group crossover,
#' the category is set to CO1 only}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Torchiano17JVLC
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Torchiano17JVLC"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Abrahao13TSE data
#'
#' Data illustrate correlations between results from individual participants in a family of five
#' cross-over experiments conducted by Abrahao et al:
#' [1] S. Abrahao, C. Gravino, E. Insfran Pelozo, G. Scanniello, and
#' G. Tortora, 'Assessing the effectiveness of sequence diagrams in the comprehension of functional
#' requirements: Results from a family of five experiments,' IEEE Transactions on Software
#' Engineering, vol. 39, no. 3, pp. 327–342, March 2013
#' The five experiments assess whether the comprehensibility of function requirements improve when
#' software models include UML sequence diagrams.
#' If you use this data set please cite:
#' [1] S. Abrahao, C. Gravino, E. Insfran Pelozo, G. Scanniello, and
#' G. Tortora, 'Assessing the effectiveness of sequence diagrams in the comprehension of
#' functional requirements: Results from a family of five experiments,' IEEE Transactions on
#' Software Engineering, vol. 39, no. 3, pp. 327–342, March 2013
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format A data frame with 224 rows and 8 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each of the five experiments in
#' the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A (DM-NODM,ECP-EPlat or MShop-Theatre ),
#'  B (NODM-DM,ECP-EPlat or  MShop-Theatre ), C(DM-NODM,EPlat-ECP or Theatre-MShop),
#'  D(NODM-DM,EPlat-ECP or Theatre-MShop)}
#' \item{System}{<fct>|Software systems used in the experiment: ECP an e-commerce platform from
#' which CDs and books can be bought, EPlat a system for the management of courses, lectures and
#' students of a university, M-Shop a system for managing sales at a music shop, Theatre a system
#' for managing bookings for a theatre.}
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment: A Dynamic Model (DM) vs No Dynamic Model (NODM)}
#' \item{Comprehension}{<dbl>|Dependent variable: The comprehension level the software engineer
#' achieved based on the F-measure }
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. For a 2 group crossover,
#' the category is set to CO1 only}
#' \item{Ability}{<fct>|Ability: An assessment of the ability of participants: Low, High, NA (not
#' available)}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Abrahao13TSE
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Abrahao13TSE"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14EASE data
#'
#' Data illustrate correlations between results from individual participants in a family of two
#' cross-over experiments conducted by Scanniello et al:
#' [1] G. Scanniello, M. Staron, H. Burden, and R. Heldal, 'On the effect of using SysML
#' requirement diagrams to comprehend requirements: results from two controlled experiments,' in
#' Proceedings of the 18th International Conference on Evaluation and Assessment in Software
#' Engineering, EASE. ACM, 2014.
#
#' The two experiments investigate whether requirements specified as SysML requirement diagrams
#' improve the comprehensibility of requirements.
#' If you use this data set please cite:
#' [1] G. Scanniello, M. Staron, H. Burden, and R. Heldal, 'On the effect of using SysML
#' requirement diagrams to comprehend requirements: results from two controlled experiments',
#' in Proceedings of the 18th International Conference on Evaluation and Assessment in Software
#' Engineering, EASE. ACM, 2014.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format A data frame with 174 rows and 9 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each experiment in the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A (RD-NORD,Automobile-ESS),
#' B (NORD-RD,ESS-Automobile), C(NORD-RD,Automobile-ESS), D(RD-NORD,ESS-Automobile).}
#' \item{System}{<fct>|Software systems used in the experiment: Automobile: A system for
#' controlling car behavior with use cases about entering the car, anti-lock breaking or operating
#' the climate control of a car. ESS (Enhanced Security System) a system designed
#' to detect potential home intruders.}
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment: RD availability of a SysML requirements diagram
#' vs No requirements diagram (NORD)}
#' \item{Time}{<dbl>|Dependent variable: The time [min] required for the comprehension task.}
#' \item{Comprehension}{<dbl>|Dependent variable: The comprehension level the software engineer
#' achieved. }
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. For 2 group crossover,
#' the category is set to CO1 only}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14EASE
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14EASE"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Ricca14TOSEM data
#'
#' Data illustrate correlations between results from individual participants in a family of three
#' of four cross-over experiments conducted by Ricca et al:
#' [1] F. Ricca, G. Scanniello, M. Torchiano, G. Reggio, and E. Astesiano,
#' 'Assessing the effect of screen mockups on the comprehension of functional requirements,'
#' ACM Transactions on Software Engineering and Methodology, vol. 24, no. 1, pp. 1:1–1:38, Oct.
#' 2014.
#' The goal of the study was to assess whether stakeholders benefit from the presence of screen
#' mock-ups in the comprehension of functional requirements represented with use cases.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format A data frame with 176 rows and 10 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each of the three experiments
#' in the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A , B , C, D}
#' \item{System}{<fct>|Software systems used in the experiment: AMICO, a system for management of
#' condominiums. EasyCoin, a system for cataloguing collections of coins.}
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment: Screen mockup available (S) vs Text only (T)}
#' \item{Time}{<dbl>|Dependent variable: The time [min] taken to perform the software engineering
#' task.}
#' \item{Comprehension}{<dbl>|Dependent variable: The comprehension level the software engineers.}
#' \item{Efficiency}{<dbl>|Dependent variable: The ratio of comprehension to time. }
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. For a 2 group crossover,
#' the category is set to CO1 only.}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Ricca14TOSEM
"KitchenhamEtAl.CorrelationsAmongParticipants.Ricca14TOSEM"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Gravino15JVLC data
#'
#' Data illustrate correlations between results from individual participants in a family of two
#' cross-over experiments conducted by Gravino et al.:
#' [1] C. Gravino, G. Scanniello, and G. Tortora, 'Source-code comprehension tasks supported by
#' UML design models: Results from a controlled experiment and a differentiated replication,'
#' Journal of Visual Languages and Computing, vol. 28, pp. 23–38, 2015.
#' The experiments assess whether the comprehension of object oriented source-code increases used
#' with UML class and sequence diagrams produced in the software design phase.
#' If you use this data set please cite:
#' [1] C. Gravino, G. Scanniello, and G. Tortora, 'Source-code comprehension tasks supported by UML
#' design models: Results from a controlled experiment and a differentiated replication,'
#' Journal of Visual Languages and Computing, vol. 28, pp. 23–38, 2015.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format A data frame with 64  rows and 9 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each of the three experiments
#' in the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A , B , C, D}
#' \item{System}{<fct>|Software systems used in the experiment: Music shop, a system for handling
#' the sales of a music shop. Theater ticket, a system for managing theatre reservations.  }
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment: Mo, design models were available, NOMo design
#' models were not available}
#' \item{Comprehension}{<dbl>|Dependent variable: The level of comprehension achieved by the
#' software engineer. }
#' \item{Time }{<dbl>|Dependent variable: The time [min] taken to complete the comprehension task. }
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. For 2 group crossover,
#' the category is set to CO1 only}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Gravino15JVLC
"KitchenhamEtAl.CorrelationsAmongParticipants.Gravino15JVLC"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14JVLC data
#'
#' Data illustrate correlations between results from individual participants in a cross-over
#' experiment conducted by Scanniello and Erra:
#' [1] G. Scanniello and U. Erra, 'Distributed modeling of use case diagrams with a method based on think-pair-square: Results from two controlled experiments', Journal of Visual Languages and
#' Computing, vol. 25, no. 4, pp. 494–517, 2014.
#' The experiment investigated whether a new method based on think-pair-square and its
#' implementation in a integrated communication/modeling environment (TPS approach) is as effective
#' as traditional face-to-face (F2F approach) for requirements elicitation. The experiment was performed in two stages using different software systems.
#' If you use this data set please cite:
#' [1] G. Scanniello and U. Erra, 'Distributed modeling of use case diagrams with a method based on
#'  think-pair-square: Results from two controlled experiments,” Journal of Visual Languages and
#'  Computing, vol. 25, no. 4, pp. 494–517, 2014.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The Importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format A data frame with 36 rows and 12 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each experiment in the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each team of four participants, unique for the specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A , B}
#' \item{System}{<fct>|Software systems used in the experiment: Library (a software system to manage
#'  books and users of a library) and FilmCollection (a software system for the selling and the
#'  rental of films in a shop) in ExperimentStage1 and Rent (a car rental software to manage
#'  cars, customers, and reservations) and ECP (an E-Commerce Platform to order CDs and books via
#'  the Internet from an on line catalogue), in ExperimentStage2. }
#' \item{Treatment}{<fct>|Experimental Treatment: TPS vs F2F.}
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2 within each stage of the
#'  experiment}
#' \item{Time}{<dbl>|Dependent variable: The total time [min] to accomplish the requirement
#' engineering task.}
#' \item{Quality}{<dbl>|Dependent variable: The quality of the requirements engineering task. }
#' \item{CrossOverID}{<fct>|Crossover category: For a single 2 group crossover experiment, the
#' value is set to CO1 for each experiment stage.}
#' \item{ExperimentPeriod}{<fct>|ExperimentPeriod: The time period across both stages of the
#'  experiment.}
#' \item{ExperimentStage}{<fct>|ExperimentStage: 1 first stage, 2 second stage.}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14JVLC
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14JVLC"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Romano18ESEM data
#'
#' Data illustrate correlations between results from individual participants in a cross-over
#' experiment conducted by Romano et al.:
#' [1] S. Romano, G. Scanniello, D. Fucci, N. Juristo, and B. Turhan, 'The effect of noise on
#' software engineers’ performance', in Proceedings of the 12th ACM/IEEE International Symposium on
#' Empirical Software Engineering and Measurement, ser. ESEM'18, 2018.
#' The experiments assess whether noise has an impact on the performance of software engineers.
#' If you use this data set please cite:
#' [1] S. Romano, G. Scanniello, D. Fucci, N. Juristo, and B. Turhan, 'The effect of noise on
#' software engineers’ performance', in Proceedings of the 12th ACM/IEEE International Symposium on
#'  Empirical Software Engineering and Measurement, ser. ESEM'18, 2018.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The Importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#' The experiment had two parts but Kitchenham et al. only use the data from the first part of the experiment.
#' @format A data frame with 194 and 10 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each part of the experiment.
#' Exp.1 identifies data from the first part of the experiment, Exp.2 identifies data from the
#' second part of the experiment.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for both
#' parts of the experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A , B }
#' \item{System}{<fct>|Software systems used in the experiment: For the first part of the
#' experiment, M-Shop (a system for managing a music shop) and Theater (a system for managing
#' theatre reservations). For the second part of the experiment: AveCalc (a system that manages as
#' electronic register and LaTazza  (a system for a drinks vending machine)  }
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment: NOISE, participants were asked to perform a
#' comprehension task in a noisy environment.
#' NORMAL, participants were asked to perform a comprehension task under normal working conditions.}
#' \item{Fc}{<dbl>|Dependent variable: the balanced F-measure which represents the trade-off
#' between  precision and recall, measured in the first part of the experiment.}
#' \item{Avg}{<dbl>|Dependent variable: The average number of fully correct answers, measured in
#' the first part of the experiment. }
#' \item{Ff}{<dbl>|Dependent variable: Effectiveness of fault correction. Measured in the second
#' part of the experiment. }
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossovers, the crossover category
#' specifies the matching pairs of sequence groups, CO1 and CO2. For 2 group crossover, the
#' category is set to CO1 only.}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Romano18ESEM
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Romano18ESEM"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Ricca10TSE data
#'
#' Data illustrate correlations between results from individual participants in a family of four
#' cross-over experiments conducted by Ricca et al.:
#' [1] F. Ricca, M. D. Penta, M. Torchiano, P. Tonella, and M. Ceccato 'How developers’ experience
#' and ability influence web application comprehension tasks supported by uml stereotypes: A series
#' of four experiments', IEEE Transactions on Software Engineering, vol. 36, no. 1, pp. 96-118,
#' 2010.
#' Although we present the full data set, only the first two experiments were used in the
#' correlation study, because many of the observations in the final two studies were unpaired.
#' The experiments assess whether participants performance comprehension tasks better when using
#' source code complemented by standard UML diagrams (UML) or by diagrams stereotyped using the
#' Conallen notation (Conallen).
#' If you use this data set please cite:
#' [1] F. Ricca, M. D. Penta, M. Torchiano, P. Tonella, and M. Ceccato 'How developers’ experience
#' and ability influence web application comprehension tasks supported by uml stereotypes: A series
#' of four experiments', IEEE Transactions on Software Engineering, vol. 36, no. 1, pp. 96—118,
#' 2010.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The Importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format A data frame with 176 rows and 10 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each of the four experiments in
#' the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A , B , C, D}
#' \item{System}{<fct>|Software systems used in the experiment: Two Java-based Web applications,
#' Claros and WfMS  }
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment: UML or Conallon}
#' \item{FMeasure}{<dbl>|Dependent variable: The comprehension level achieved by the participant.}
#' \item{Time}{<dbl>|Dependent variable: The time [min] to complete the experimental task}
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. For 2 group crossover,
#' the category is set to CO1 only}
#' \item{Ability}{<fct>| h: High l: Low, NA: Not available}
#' \item{Experience}{<fct>| G: Master students, U: undergraduates, P: researchers}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Ricca10TSE
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Ricca10TSE"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Reggio15SSM data
#'
#' Data illustrate correlations between results from individual participants in a family of two
#' cross-over experiments conducted by Reggio et al:
#' [1] G. Reggio, F. Ricca, G. Scanniello, F. D. Cerbo, and G. Dodero,'On the comprehension of
#' workflows modeled with a precise style: results from a family of controlled experiments'.
#' Software and Systems Modeling, vol. 14, pp. 1481–1504, 2015.
#' The experiments assess whether the level of formality/precision in workflow model influences
#' comprehension.
#' If you use this data set please cite:
#' [1] G. Reggio, F. Ricca, G. Scanniello, F. D. Cerbo, and G. Dodero, 'On the comprehension of
#' workflows modeled with a precise style: results from a family of controlled experiments'.
#' Software and Systems Modeling, vol. 14, pp. 1481–1504, 2015.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'The Importance
#' of the Correlation between Results from Individual Participants in Crossover Experiments'
#' (to be submitted as of 2020).
#'
#' @format A data frame with 78 rows and 9 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each of the three experiments
#' in the data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A , B , C, D}
#' \item{System}{<fct>|Software systems used in the experiment: PO, a system to process orders for
#' an online shop. DM, a system to manage an online document review process.}
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment: }
#' \item{Comprehension}{<dbl>|Dependent variable: The comprehension level obtained by each
#' participant.}
#' \item{Time}{<dbl>|Dependent variable: The time [min] taken by each participant to complete the
#' comprehension task.}
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. For a 2 group crossover,
#' the category is set to CO1 only}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Reggio15SSM
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Reggio15SSM"


#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello17TOSEM data
#'
#' Data illustrate correlations between results from individual participants in a family of four
#' cross-over experiments conducted by Scanniello et al.:
#' [1] G. Scanniello, M. Risi, P. Tramontana, and S. Romano, 'Fixing faults in C and Java source
#' code: Abbreviated vs. full-word identifier names', ACM Transactions on Software Engineering
#' Methodology, vol. 26, no. 2, 2017.
#' The experiments assess whether whether the use of abbreviated identifier names (ABBR), impacts
#' the effectiveness of fault fixing in C and Java source code in comparison with full-word
#' identifier names (FULL).
#' If you use this data set please cite:
#' [1] G. Scanniello, M. Risi, P. Tramontana, and S. Romano, “Fixing faults in C and Java source
#' code: Abbreviated vs. full-word identifier names', ACM Transactions on Software Engineering
#' Methodology, vol. 26, no. 2, 2017.
#' [2] Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and Carmine Gravino, 'On the
#' Importance of the Correlation between Results from Individual Participants in Crossover
#' Experiments' (to be submitted as of 2020).
#'
#' @format A data frame with 200 rows and 17 variables:
#' \describe{
#' \item{ExperimentID}{<fct>|ExperimentID: A unique identifier for each of the experiments in the
#' data set.}
#' \item{ParticipantID}{<fct>|Participant ID: An identifier for each participant, unique for a
#' specific experiment.}
#' \item{SequenceGroup}{<fct>|Experimental Sequence Group: A , B , C, D}
#' \item{System}{<fct>|Software systems used in the experiments: The Unibas experiment used Agenda
#' (a system for tracking personal contacts) and Gas-Station (a system for managing a petrol
#' station). The UniNa experiment used Financial (a system which is a command line option price
#' calculator) and Hotel-Reservation. The POLINA and PROF experiments used AveCalc (a system that
#' manages as electronic register and LaTazza  (a system for a drinks vending machine).}
#' \item{Period}{<fct>|Time period of the cross-over experiment: 1 or 2}
#' \item{Treatment}{<fct>|Experimental Treatment:ABBR, abbreviated names. FULL, full names}
#' \item{Time}{<dbl>|Dependent variable: The time each participant spent performing the SE task.}
#' \item{FMeasure}{<dbl>|Dependent variable: The effectiveness of the participants taking into
#' account correctness and completeness of the fault fixing tasks}
#' \item{Efficiency}{<dbl>|Dependent variable: The ratio of effectiveness to time. }
#' \item{CrossOverID}{<fct>|CrossOver category: For 4 group crossover designs, the crossover
#' category specifies the matching pairs of sequence groups, CO1 and CO2. For 2 group crossover,
#' the category is set to CO1 only.}
#' \item{Language}{<fct>|Java or C. The language was the same for all participants in a specific
#' experiment. POLINA and PROF used Java, UNIBAS and UNINA used C.}
#' \item{Ident}{<dbl>|Dependent variable: The number of faults identified.}
#' \item{Fixed}{<dbl>|Dependent variable: The number of faults identified.}
#' \item{WrongIdent}{<dbl>|Dependent variable: The number of faults incorrectly identified}
#' \item{WronglyFixed}{<dbl>|Dependent variable: The number of faults incorrectly fixed.}
#' \item{precision}{<dbl>|Dependent variable: The ratio of number of faults correctly fixed to the
#' number of faults correctly identified.}
#' \item{recall}{<dbl>|Dependent variable: The ratio of number of faults correctly fixed to the
#' total number of fault.}
#' }
#' @examples
#' KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello17TOSEM
#'
"KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello17TOSEM"
