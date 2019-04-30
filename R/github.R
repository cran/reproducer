#' @title searchForIndustryRelevantGitHubProjects
#' @description Function searches for industry relevant software projects available from GitHub. The function was used to deliver data set of software projects in an NCBiR project. More details are desribed in a report: Lech Madeyski, “Training data preparation method,” tech. rep., code quest (research project NCBiR POIR.01.01.01-00-0792/16), 2019, as well as a paper: Tomasz Lewowski and Lech Madeyski, "Creating evolving project data sets in software engineering", 2019.
#' If you use this function or the returned data set than please cite:
#' Tomasz Lewowski and Lech Madeyski, "Creating evolving project data sets in software engineering", 2019
#' @author Lech Madeyski and Tomasz Lewowski
#' @export searchForIndustryRelevantGitHubProjects
#' @param myToken A private token used to access GitHub
#' @return selected GitHub projects
#' @examples
#'   #to run this function you need to use your own token as a parameter of the function
#'   #calculateSmallSampleSizeAdjustment("...") #use your own token as a parameter of the function
#' @importFrom dplyr %>%
searchForIndustryRelevantGitHubProjects = function(myToken) {
    token <- myToken

    language <- "Java"
    minForks <- 1
    minStars <- 1

    pushedAt <- "2018-09-01..2019-03-25"
    repoCount <- 50
    languageCount <- 1
    afterClause <- ''
    listOfRepos <- NULL
    repoCountTotal <- 100
    analyzedRepos <- repoCountTotal / repoCount
    hasNextPage <- TRUE

    companies <-
    c(
    "microsoft",
    "google",
    "redhat-developer",
    "IBM",
    "Intel",
    "amzn",
    "SAP",
    "thoughtworks",
    "alibaba",
    "github",
    "facebook",
    "Tencent",
    "Pivotal",
    "epam",
    "baidu",
    "mozilla",
    "oracle",
    "Unity-Technologies",
    "uber",
    "yandex",
    "Shopify",
    "linkedin",
    "SUSE",
    "Esri",
    "Apple",
    "salesforce",
    "vmware",
    "adobe",
    "andela",
    "cisco"
    )
    # Who really contributes to open source https://www.infoworld.com/article/3253948/open-source-tools/who-really-contributes-to-open-source.html


    AmazonGitHubAccountsFrom_amzn.github.io <-
    c(
    "alexa",
    "aws",
    "awsdocs",
    "awslabs",
    "aws-quickstart",
    "aws-samples",
    "blox",
    "boto",
    "carbonado",
    "c9",
    "ajaxorg",
    "cloud9ide",
    "goodreads",
    "gluon-api",
    "IvonaSoftware",
    "twitchtv",
    "twitchtv",
    "TwitchScience",
    "justintv",
    "Zappos",
    "amazon-archives"
    )

    #Pivotal also keep projects under project name, not company (vide https://pivotal.io/open-source )
    #Pivotal is involved in or supports some Apache projects...
    OtherPivotalGitHubAccounts <-
    c(
    "pivotalsoftware",
    "Pivotal-Open-Source-Hub",
    "kubernetes",
    "spring-projects",
    "cloudfoundry",
    "reactor",
    "openservicebrokerapi",
    "greenplum-db",
    #Apache and Exlipse projects
    "rabbitmq",
    "jasmine",
    "robolectric",
    "pivotal-sprout"
    )


    Apache <- "apache"

    Eclipse <- "eclipse"

    companies <-
    c(
    companies,
    AmazonGitHubAccountsFrom_amzn.github.io,
    OtherPivotalGitHubAccounts,
    Apache,
    Eclipse
    )

    for (company in companies)
    {
        currentPage <- 0
        afterClause <- ''
        hasNextPage <- TRUE
        searchQuery <-
        GetoptLong::qq(
        "\\\"user:@{company} language:@{language}\\\"")
        print(searchQuery)

        while (hasNextPage) {
            query <- gsub("\n", " ", GetoptLong::qq(
            '{"query": "query {
            viewer {
            login
            }
            rateLimit {
            limit
            cost
            remaining
            resetAt
            }
            search(query:@{searchQuery}, type: REPOSITORY, first: @{repoCount}@{afterClause}){
            repositoryCount
            pageInfo {
            endCursor
            hasNextPage
            }
            nodes {
            ... on Repository {
            id
            nameWithOwner
            createdAt
            updatedAt
            pushedAt
            diskUsage
            forkCount
            isArchived
            isFork
            isLocked
            isMirror
            isPrivate
            sshUrl
            licenseInfo {
            name
            }
            defaultBranchRef {
            target {
            ... on Commit {
            oid
            history {
            totalCount
            }
            }
            }
            }
            stargazers {
            totalCount
            }
            collaborators {
            totalCount
            }
            watchers {
            totalCount
            }
            languages(first: @{languageCount}, orderBy:{direction:DESC, field: SIZE}) {
            totalSize
            edges {
            size
            node {
            name
            }
            }
            }
            }
            }
            }
        }","variables": {}}
            '
            ))

            #  cli$load_schema()

            result <- httr::POST(
              url="https://api.github.com/graphql",
              httr::add_headers(
                  Authorization = paste0("Bearer ", token),
                  Accept = "application/vnd.github.machine-man-preview+json"
                  ),
              body=query
            )
            #qry <- ghql::Query$new()
            #qry$query('myquery', query)

            jsonResult <- httr::content(result)
            repos <- jsonlite::fromJSON(jsonlite::toJSON(jsonResult), flatten = TRUE)

            rows <- repos$data$search$nodes
            rows["searchQuery"] <- searchQuery

            listOfRepos <- dplyr::bind_rows(listOfRepos, rows)

            repoCountTotal <- repos$data$search$repositoryCount

            hasNextPage <- repos$data$search$pageInfo$hasNextPage
            print(
            GetoptLong::qq(
            "Page: @{currentPage + 1} ends at cursor: @{repos$data$search$pageInfo$endCursor}. Has next page? @{hasNextPage}"
            )
            )
            print(
            GetoptLong::qq(
            "In total there are: @{repos$data$search$repositoryCount} repositories. This query contained: @{repoCount} starting from @{currentPage * repoCount}"
            )
            )
            print(
            GetoptLong::qq(
            "Query cost: @{repos$data$rateLimit$cost}, remaining: @{repos$data$rateLimit$remaining}/@{repos$data$rateLimit$limit} until: @{repos$data$rateLimit$resetAt} for @{repos$data$viewer$login}"
            )
            )
            afterClause <-
            GetoptLong::qq(", after: \\\"@{repos$data$search$pageInfo$endCursor}\\\"")
            currentPage <- currentPage + 1
        }


        print(
        GetoptLong::qq(
        "Finished fetching after @{currentPage} pages. Downloaded metadata for @{currentPage * repoCount} repositories. Day: @{pushedAt}"
        )
        )

    }

    listOfRepos <- listOfRepos %>% tidyr::drop_na("id")

    allGithubProjects <- listOfRepos

    allGithubProjects <-
    tibble::rowid_to_column(allGithubProjects, "rowID")

    allGithubProjects <-
    allGithubProjects %>% dplyr::rename(commitSHA = "defaultBranchRef.target.oid")
    allGithubProjects <-
    allGithubProjects %>% dplyr::rename(sshUrlOfRepository = "sshUrl")

    # split column "languages.edges" (including, e.g., "1104419, Java") into two columns "Java.byte.count" (inclunding, e.g., "1104419") and "Language" (e.g., "Java")
    allGithubProjects <-
    tidyr::separate(
    allGithubProjects,
    col = "languages.edges",
    into = c("Java.byte.count", "Language"),
    convert = TRUE,
    sep = "\\,"
    )
    allGithubProjects$Java.byte.count <-
    readr::parse_number(allGithubProjects$Java.byte.count)
    #https://stackoverflow.com/questions/29508943/r-regular-expression-isolate-a-string-between-quotes
    allGithubProjects$Language <-
    sub('[^\"]+\"([^\"]+).*', '\\1', "Language")

    myselection <- function(projects){
        stringr::str_sub(projects$createdAt, 1, 10) < "2018-01-01" &
            stringr::str_sub(projects$pushedAt, 1, 10) > "2018-09-01" &
            projects$isArchived == FALSE &
            projects$watchers.totalCount > stats::quantile(unlist(projects$watchers.totalCount), 1 /
                                                               4) &
            projects$stargazers.totalCount > stats::quantile(unlist(projects$stargazers.totalCount), 1 /
                                                                 4) &
            projects$Java.byte.count > stats::quantile(unlist(projects$Java.byte.count), na.rm = TRUE, 1 /
                                                           4) &
            projects$defaultBranchRef.target.history.totalCount > stats::quantile(
                unlist(projects$defaultBranchRef.target.history.totalCount),
                1 / 4
            ) &
            projects$forkCount > stats::quantile(unlist(projects$forkCount), 1 / 4)

    }

    selectedGithubProjects <-
    subset(
    allGithubProjects,
    myselection(allGithubProjects)
    )


    sshPaths2SelectedProjects <-
    subset(
    allGithubProjects,
    myselection(allGithubProjects),
    select = c("rowID", "sshUrlOfRepository", "commitSHA")
    )

    sshPaths2SelectedProjects <-
    tibble::rowid_to_column(sshPaths2SelectedProjects, "rowNr")


    sshPaths2AllProjects <-
    subset(allGithubProjects,
    select = c("rowID", "sshUrlOfRepository", "commitSHA"))

    repoOwners <-
    dplyr::count(selectedGithubProjects,
    stringr::word("nameWithOwner", 1, sep = "\\/"))

    openxlsx::write.xlsx(
    as.data.frame(allGithubProjects),
    file = paste("githubAllProjects", ".xlsx", sep = ""),
    colNames = TRUE,
    borders = "columns"
    )
    openxlsx::write.xlsx(
    as.data.frame(selectedGithubProjects),
    file = paste("githubSelectedProjects", ".xlsx", sep = ""),
    colNames = TRUE,
    borders = "columns"
    )

    openxlsx::write.xlsx(
    as.data.frame(sshPaths2SelectedProjects),
    file = paste("githubSshPaths2SelectedProjects", ".xlsx", sep = ""),
    colNames = TRUE,
    borders = "columns"
    )

    return(selectedGithubProjects)
}


