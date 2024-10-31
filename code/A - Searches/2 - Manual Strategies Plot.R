# CONTEXT 
# 
# Following manual review, the protocols were summarized according to:
# 1. Explanation 
#    Whether damaged filtering protocol was explained 
#
# 2. Strategy used 
#    If explained, what strategy was used. First strategy was matched to of the following 11  
#    but if no appropriate match was found the method was added to other 
# 
#    Tools: ddqc, DropletQC, ensembleKQC, miQC, scater, valiDrops
#    Manaul (outlier detection of any kind): 
#     i. All
#        Any outlier based filtering based on 5 or more criteria such as feature and UMI counts,
#        mitochondrial, ribosomal and MALAT1 expression/ percentages.
#    ii. Mito -ribo 
#        Any outlier based filtering focusing mitochondrial and ribosomal expression/ 
#        percentages, may include feature and UMI counts.
#   iii. Mito 
#        Any outlier based filtering focusing mitochondrial expression/ 
#        percentages, may include feature and UMI counts.
#    iv. Mito -isolated
#        Any outlier based filtering focusing on only mitochondrial expression/ 
#        percentages and nothing else
#     v. MALAT1
#        Any outlier based filtering focusing on only MALAT1 expression/ 
#        percentages and nothing else

reviewed_methods <- read.csv("/home/alicen/Projects/ReviewArticle/article_search/PubMed_search_papers_reviewed.csv")

