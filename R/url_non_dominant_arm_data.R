# url_non_dominant_arm_data
fmsURL<-"http://www.stat-gen.org/book.e1/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)

# see page 21 of Applied Statistical Genetics with R. Andrea S. Foulkes
