msBayes: Pipeline for testing comparative Phylogeographic histories
using hierarchical ABC

Pronounced "em es bayeszzz"

* [Installation Instructions](https://docs.google.com/document/d/1enMQaogxOs0RppAmE8KcGU3nNjzotuiycAl6I1s0KYg/edit)
* [User Guide](https://docs.google.com/document/d/15heQlz60cGe6GWKcXqf1AIYMBZ6p2sKuGct-EoEHiNU/edit)

## Dockerfile

In order to allow users to run this tools in an easy way was added to the repository a Dockerfile that allows you to run this tool in any host with docker installed.

### Prerequisits
- Linux Machine
- Docker service installed
- git installed

### Running ABC with docker

1. *Clone Repo:* first you need to clone this repo.
2. *Build Docker Image:* for build docker image you need to get into repo directory, then run the following command.
```
docker build -t msbayes:latest .
```
3. Run the first step of the procedure (for the example in the example folder it can't be run due missing .im files)

```
docker run --rm -v /path/to/myworkspace:/workspace --name msbayes msbayes:latest convertIM.pl infile.list
```
**/path/to/myworkspace:** Path to your working directory where your have your files for work with msbayes
**infile.list:** Text plane list of your *.im files relative to your working directory

4. Run second step.

```
docker run --rm -v /home/ubuntu/workspace:/workspace --name msbayes msbayes:latest obsSumStats.pl -T obsSS.table batch.masterIn.fromIM > obsSS.txt
```

5. Run step three

```
docker run --rm -v /home/ubuntu/workspace:/workspace --name msbayes msbayes:latest msbayes.pl -c batch.masterIn.fromIM -s 7 -r 5000 -o priorfile
```
6. Run step four

```
docker run --rm -v /home/ubuntu/workspace:/workspace --name msbayes msbayes:latest acceptRej.pl -p outfig.pdf obsSS.txt priorfile > modeEstimatesOut.txt
```

7. Enjoy the results :)

