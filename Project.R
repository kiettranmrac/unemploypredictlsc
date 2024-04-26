###### Sample data
# E is the matrix of the data from 2009-2019 in months, C is the training data from 2009-2014 and D is the "soon to be" generated data
E = c (7.8, 8.3, 8.7, 9.0,	9.4,	9.5,	9.5,	9.6,	9.8,	10.0, 9.9, 9.9,
9.8, 9.8,	9.9,	9.9,	9.6,	9.4,	9.4,	9.5,	9.5,	9.4,	9.8,	9.3,
9.1, 9.0,	9.0,	9.1,	9.0,	9.1,	9.0,	9.0,	9.0,	8.8,	8.6,	8.5,
8.3, 8.3,	8.2,	8.2,	8.2,	8.2,	8.2,	8.1,	7.8,	7.8,	7.7,	7.9,
8.0, 7.7,	7.5,	7.6,	7.5,	7.5,	7.3,	7.2,	7.2,	7.2,	6.9,	6.7,
6.6, 6.7,	6.7,	6.2,	6.3,	6.1,	6.2,	6.1,	5.9,	5.7,	5.8,	5.6,
5.7, 5.5,	5.4,	5.4,	5.6,	5.3,	5.2,	5.1,	5.0,	5.0,	5.1,	5.0,
4.8, 4.9,	5.0,	5.1,	4.8,	4.9,	4.8,	4.9,	5.0,	4.9,	4.7,	4.7,
4.7, 4.6,	4.4,	4.4,	4.4,	4.3,	4.3,	4.4,	4.3,	4.2,	4.2,	4.1,
4.0, 4.1,	4.0,	4.0,	3.8,	4.0,	3.8,	3.8,	3.7,	3.8,	3.8,	3.9,
4.0,3.8, 3.8, 3.7, 3.6, 3.6, 3.7, 3.6, 3.5, 3.6, 3.6, 3.6)
dim ( E ) = c (72, 1)
E
C = c (7.8, 8.3, 8.7, 9.0,	9.4,	9.5,	9.5,	9.6,	9.8,	10.0, 9.9, 9.9,
9.8, 9.8,	9.9,	9.9,	9.6,	9.4,	9.4,	9.5,	9.5,	9.4,	9.8,	9.3,
9.1, 9.0,	9.0,	9.1,	9.0,	9.1,	9.0,	9.0,	9.0,	8.8,	8.6,	8.5,
8.3, 8.3,	8.2,	8.2,	8.2,	8.2,	8.2,	8.1,	7.8,	7.8,	7.7,	7.9,
8.0, 7.7,	7.5,	7.6,	7.5,	7.5,	7.3,	7.2,	7.2,	7.2,	6.9,	6.7,
6.6, 6.7,	6.7,	6.2,	6.3,	6.1,	6.2,	6.1,	5.9,	5.7,	5.8,	5.6)
dim ( C ) = c (72, 1)
C
D = c (5.7, 5.5,	5.4,	5.4,	5.6,	5.3,	5.2,	5.1,	5.0,	5.0,	5.1,	5.0,
4.8, 4.9,	5.0,	5.1,	4.8,	4.9,	4.8,	4.9,	5.0,	4.9,	4.7,	4.7,
4.7, 4.6,	4.4,	4.4,	4.4,	4.3,	4.3,	4.4,	4.3,	4.2,	4.2,	4.1,
4.0, 4.1,	4.0,	4.0,	3.8,	4.0,	3.8,	3.8,	3.7,	3.8,	3.8,	3.9,
4.0,3.8, 3.8, 3.7, 3.6, 3.6, 3.7, 3.6, 3.5, 3.6, 3.6, 3.6)
dim ( D ) = c (60, 1)
D

###### Function to calculate the "best fit" line
UF = function (C, D, ct ) {
	# Given a C matrix, generate the matrix X with items as
	# u_1,u_2,u_3...
	# u_2,u_3,u_4...
	# ...
	# ...u_n-2,u_n-1,u_n
	X = array (dim = c (dim(C)[1] - ct, ct))
	i=1
	while ( i <= dim(C)[1] - ct) {
		j = 1
		while (j <= ct){
			X[i,j] = C[i + j - 1,1]
			j= j + 1
		}
		i = i + 1
	}
	# X # uncomment to see the corresponding value
	
	# Calculate the y vector
	y = array (dim = c (dim(X)[1], 1))
	i=1
	while ( i <= dim(y)[1]) {
		y[i,1] = C[i + ct,1]
		i = i + 1
	}
	# y  # uncomment to see the corresponding value

	# Calculate the number of "offsprings" (the top values) for the Leslie models
	XT = t(X) # X^T
	# XT # uncomment to see the corresponding value
	
	XTX  = XT%*%X # X^T * X
	#XTX # uncomment to see the corresponding value
	
	XTXI = solve(XTX) # inverse of X^T * X
	# XTXI # uncomment to see the corresponding value
	
	XTy = XT%*%y # X^T * y
	# XTy # uncomment to see the corresponding value
	
	b = XTXI%*%XTy # (inverse of X^T*X) * X^T * y
	# b # uncomment to see the corresponding value
	
	# Generate the matrix for the Leslie model
	S = array (dim = c (dim(b)[1], dim(b)[1]))
		# The top values
		i=1
		while ( i <= dim(b)[1]) {
			S[1, dim(S)[1]- i+1] = b[i, 1]
			i = i + 1
		}
		
		# The other values (assumed to be 1s for the proporional of survival
		# since we want the "full" effect of the past data")
		i=1
		while ( i <= dim(b)[1] - 1) {
			S[i+1,i] = 1
			i = i + 1
		}
		S[is.na(S)] <- 0
	# S # uncomment to see the corresponding value
	eigen(S)$values # uncomment to see the corresponding value
		
	# Calculate the value for the other years
		# First, calculate the initial vector
		init = array (dim = c(ct, 1))
		i=1
		while ( i <= ct) {
			init[ct + 1 - i, 1] = D[i,1]
			i = i + 1
		}
		# init # uncomment to see the corresponding value
	
		# Second, calculate the other years and save it in a vector
		res = array (dim = c(dim(D)[1], 1))
		i = 1
		while ( i <= dim(res)[1]) {
			init = S%*%init
			res[i, 1] = init[1]
			i = i + 1
		}
	res
}

###### Calculate the errors for the usage of past datas
error = array (dim = c (dim(D)[1],1))
ct = 1 # Order of the model
while ( ct <= dim(C)[1]/2) {
		res = UF(C, D, ct)
		res
		D
		error[ct,1] = sum((D-res)^2)/(dim(D)[1]-ct)
		ct = ct + 1
}
error

############################################# GRAPH AND COMPARE

###### Graph the error vector
TIME = c (1: dim((error))[1])
plot ( TIME , error, ylim = c (0 , max (10)) , xlab = "Dimension" , ylab = "Error" , pch = 19)
title ( main = "Error")
lines ( TIME , error) 
points ( TIME , error , pch = 19 , col = " red ")

###### The graph is a comparison between actual data in red and generated data in blue
TIME = c (1: dim(res)[1])
plot ( TIME , res, ylim = c (0 , max ( D)) , xlab = "Time" , ylab = "Unemployment Rate" , pch = 19)
title ( main = "Unemployment Rate, red: actual data, blue: generated data")
lines ( TIME , res)
points ( TIME , res , pch = 19 , col = " blue ")
lines ( TIME , D)
points ( TIME , D, pch = 19 , col = " red ")

