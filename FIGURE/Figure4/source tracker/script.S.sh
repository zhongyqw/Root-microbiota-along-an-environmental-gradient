#!/bin/sh

#......
#echo "" >> ~/.bash_profile
#echo "export SOURCETRACKER_PATH=/home/zhongyangqw/sourcetracker-master/project_rootmb2/" >> ~/.bash_profile 
#source ~/.bash_profile

# get help:
# Rscript sourcetracker_for_qiime.r-h 
#
# run sink predictions using QIIME taxon abundance file:
# Rscript sourcetracker_for_qiime.r-t taxa.txt -m map.txt 
#
# run leave-one-out source-sample predictions using QIIME taxon abundance file:
# Rscript sourcetracker_for_qiime.r-t taxa.txt -m map.txt -s 
#
# run sink predictions using QIIME OTU table:

#mkdir R.sink.result
#Rscript sourcetracker_for_qiime.r -i /home/zhongyangqw/sourcetracker-master/project_rootmb2/data/otus.txt -m /home/zhongyangqw/sourcetracker-master/project_rootmb2/data/metadata.R.sink.txt -o /home/zhongyangqw/sourcetracker-master/project_rootmb2/R.sink.result/

#mkdir RB.sink.result
#Rscript sourcetracker_for_qiime.r -i /home/zhongyangqw/sourcetracker-master/project_rootmb2/data/otus.txt -m /home/zhongyangqw/sourcetracker-master/project_rootmb2/data/metadata.RB.sink.txt -o /home/zhongyangqw/sourcetracker-master/project_rootmb2/RB.sink.result/

#mkdir RJ.sink.result
#Rscript sourcetracker_for_qiime.r -i /home/zhongyangqw/sourcetracker-master/project_rootmb2/data/otus.txt -m /home/zhongyangqw/sourcetracker-master/project_rootmb2/data/metadata.RJ.sink.txt -o /home/zhongyangqw/sourcetracker-master/project_rootmb2/RJ.sink.result/

mkdir S.sink.result
Rscript sourcetracker_for_qiime.r -i /home/zhongyangqw/sourcetracker-master/project_rootmb2/data/otus.txt -m /home/zhongyangqw/sourcetracker-master/project_rootmb2/data/metadata.S.sink.txt -o /home/zhongyangqw/sourcetracker-master/project_rootmb2/S.sink.result/
