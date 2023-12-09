##### Similarity Score Functions #####
# Note that the input for all similarity measures are two 1-d arrays of the same length. 
# These 1-d arrays must be normalized to sum to 1 for the Renyi and Tsallis Entropy Similarity Measures.


#Cosine Similarity Measure
S_cos = function(X,Y){
  (X %*% Y) / (sqrt(sum(X^2)) * sqrt(sum(Y^2)))
}


#Compute the Shannon entropy of a probability distribution
ent_shannon = function(vec){
  return(-sum(vec*log(vec)))
}

#Compute the Renyi entropy of a probability distribution
ent_renyi = function(vec){
  return(log(sum(vec^q)) / (1-q))
}

#Compute the Tsallis entropy of a probability distribution
ent_tsallis = function(vec){
  return(sum((vec^q)-1) / (1-q))
}


#Shannon Entropy Similarity Measure
#This similarity function was presented by: 
#Li, Y.; Kind, T.; Folz, J.; Vaniya, A.; Mehta, S. S.; Fiehn, O.
#Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. 
#Nature Methods 2021, 18, 1524–1531
S_shannon = function(X,Y){
  ent_x = ent_shannon(X)
  ent_y = ent_shannon(Y)
  ent_xy = ent_shannon((X+Y)/2)
  return(1 - (2*ent_xy - ent_x - ent_y) / log(4))
}


#Renyi Entropy Similarity Measure
#This is a novel similarity measure which generalizes the Shannon Entropy Similarity Measure
#Note that the Renyi Similarity Measure approaches the Shannon Entropy Similiarity Measure as q approaches 1
#Note that X and Y must be normalized to sum to 1
S_renyi = function(X,Y,q){
  ent_x = ent_renyi(X)
  ent_y = ent_renyi(Y)
  ent_xy = ent_renyi((X+Y)/2)
  N = (1/(1-q)) * (2*log(sum((X/2)^q) + sum((Y/2)^q)) - log(sum(X^q)) - log(sum(Y^q)))
  return(1 - (2*ent_xy - ent_x - ent_y) / N)
}


#Tsallis Entropy Similarity Measure
#This is a novel similarity measure which generalizes the Shannon Entropy Similarity Measure
#Note that the Tsallis Similarity Measure approaches the Shannon Entropy Similiarity Measure as q approaches 1
#Note that X and Y must be normalized to sum to 1
S_tsallis = function(X,Y,q){
  ent_x = ent_tsallis(X)
  ent_y = ent_tsallis(Y)
  ent_xy = ent_tsallis((X+Y)/2)
  N = sum(2*(X/2)^q + 2*(Y/2)^q - X^q - Y^q) / (1-q)
  return(1 - (2*ent_xy - ent_x - ent_y) / N)
}







