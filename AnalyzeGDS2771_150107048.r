# 1. Install packages to read the NCBI's GEO microarray SOFT files in R
# 1.Ref. http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/

# 1.1. Uncomment only once to install stuff

#source("https://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("Affyhgu133aExpr")


# 1.2. Use packages # Comment to save time after first run of the program in an R session

library(Biobase)
library(GEOquery)

# Add other libraries that you might need below this line



# 2. Read data and convert to dataframe. Comment to save time after first run of the program in an R session
# 2.1. Once download data from ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS2nnn/GDS2771/soft/				
# 2.Ref.1. About data: http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS2771
# 2.Ref.2. Study that uses that data http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694402/pdf/nihms471724.pdf
# 2.Warning. Note that do not use FULL SOFT, only SOFT, as mentioned in the link above. 2.2.R. http://stackoverflow.com/questions/20174284/error-in-gzfilefname-open-rt-invalid-description-argument

gds2771 <- getGEO(filename='C:/Users/User/Downloads/GDS2771.soft.gz') # Make sure path is correct as per your working folder. Could be './GDS2771.soft.gz'
eset2771 <- GDS2eSet(gds2771) # See http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/

# 2.2. View data (optional; can be commented). See http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/
eset2771 # View some meta data
featureNames(eset2771)[1:10] # View first feature names
sampleNames(eset2771) # View patient IDs. Should be 192
pData(eset2771)$disease.state #View disease state of each patient. Should be 192
	
# 2.3. Convert to data frame by concatenating disease.state with data, using first row as column names, and deleting first row
data2771 <- cbind2(c('disease.state',pData(eset2771)$disease.state),t(Table(gds2771)[,2:194]))
colnames(data2771) = data2771[1, ] # the first row will be the header
data2771 = data2771[-1, ] 

# 2.4. View data frame (optional; can be commented)
View(data2771)



# WRITE YOUR CODE BELOW THIS LINE


x=data2771[,-1]

#We will check if the output is a matrix or not
#is.matrix(x)
#Now the output is matrix so we will check its type
str(x)
#The matrix is a character type matrix
#We will convert character matrix to numerical matrix
storage.mode(x) <- "numeric"
#Here we are removing all the coulmns containing NA
x <- x[,colSums(is.na(x)) != nrow(x)]
#Now we will check if the matrix is numeric or not
str(x)
#is.numeric(x)  #returns true
y=data2771[,1] #taking the factors in y
#y <- as.numeric(y) #coverting the factors to numeric
is.matrix(y)  #returns false
y=as.matrix(y,nrow(data2771),1)# converting y into a 192 * 1 numeric matrix


y=ifelse(y=="3","2",y) #We will assume that the suspect cancer as no cancer


#install.packages("glmnet", repos = "http://cran.us.r-project.org")    #uncomment to install 


library("glmnet") #Loading glmnet library, Comment it after running it first time
#?glmnet   #uncomment to know use of glmnet
s=sample(192,40)#here we are using 40 samples for Cross validation

#apply the same thing for different alpha values. Here we are using nfold cv as 7
fit=cv.glmnet(x[-s,],y[-s,],family="binomial",type.measure="class",alpha=0,nfolds=7)

#the following funtion is calculating error
z=predict(fit,x[s,],s=fit$lambda.1se)#generally lambda.1se gives the best results
#summary(z)
#summary(fit)
#plot(fit)
error_analyse_cv <- function(val_thresold ,s , z) #val_thresold is for sending the value of expected threshold, s is the sample data to be sent as cv. It is same as previous s. z the prdicted value of cross validation
{# here we are calculting the percent accuracy for neighbour hood of theta. This will help us in analysing the data.To view the data uncomment the plot command
#you can assume the initial val_thresold as 0 for good results. The best val_thresold will lie on this neighbourhood most likely 
percent_accuracy=0
sm3=matrix(0,21,1)#storing percent error for different x[j]
val3=matrix(0,21,1)#storing different values of val_thresold

for(j in 1:21)
	{	
	
	pred=ifelse(z >= val_thresold+((j-11)/5)*.1,"2","1")#calculating the expected predicted matrix for different val . It will calculate (-.2,+.2) neighbour of 

	K=y[s ] #here we are storing the value of total number of correct output
	sm=0  #here sm calculates totalnumber of correct prediction for specific val_thresold
	for (i in seq_len(40))
		{
		if (pred[i] == K[i])
			{
			sm=sm+1
			}
		}	
	sm3[j]=sm/40
	val3[j]= val_thresold+((j-11)/5)
	}
plot(val3,sm3)	#run the plot few times to get the expected val from the graph produced 
}


#use this funtion to analyse what should be the best thresold to estimate results
efficiency_cal<-function(eff_val,s, z)#here eff_val is the thresold,here s is sample inputs for CV
{
	pred=ifelse(z >= eff_val,"2","1")#calculating the expected predicted matrix for different val . It will calculate efficiency at specific alpha and thresold
	K=y[s ] #here we are storing the value of total number of correct output
	sm=0  #here sm calculates totalnumber of correct prediction for specific val_thresold
	for (i in seq_len(40))
		{
		if (pred[i] == K[i])
			{
			sm=sm+1
			}
		}	
	return(sm*100/40)#returns percent accuracy
	}


	
	
avg_percent<-function(alpha1,eff_val)#this funtion calculates average percentage for a given data values
{
percent_acc=matrix(0,10,1)
for (i in 1:10)
		{
		s=sample(192,10);
		fit=cv.glmnet(x[-s,],y[-s,],family="binomial",type.measure="class",alpha=alpha1,nfolds=7)#I am using 7 fold CV
		z=predict(fit,x[s,],s=fit$lambda.1se)#generally lambda.1se gives the best results
		percent_acc[i]=efficiency_cal(eff_val,s,z)
		}
return(c(max(percent_acc),sum(percent_acc)/10))#returns max and average of percent
}


#now source this file and use different alphas and compute different values to thresold 
#you have to use s=sample(192,40) many times to many values more real.	


