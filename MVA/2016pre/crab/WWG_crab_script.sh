#this is not mean to be run locally
#
echo Check if TTY
if [ "`tty`" != "not a tty" ]; then
    echo "YOU SHOULD NOT RUN THIS IN INTERACTIVE, IT DELETES YOUR LOCAL FILES"
else

    echo "ENV..................................."
    env 
    echo "VOMS"
    voms-proxy-info -all
    echo "CMSSW BASE, python path, pwd"
    echo $CMSSW_BASE 
    echo $PYTHON_PATH
    echo $PWD 
    rm -rf $CMSSW_BASE/lib/
    rm -rf $CMSSW_BASE/src/
    rm -rf $CMSSW_BASE/module/
    rm -rf $CMSSW_BASE/python/
    mv lib $CMSSW_BASE/lib
    mv src $CMSSW_BASE/src
    mv module $CMSSW_BASE/module
    mv python $CMSSW_BASE/python

    isdata=""
    year=""
    era=""
    for i in "$@"
    do
        case $i in
            isdata=*)
            isdata="${i#*=}"
            ;;
        esac
        case $i in
            year=*)
            year="${i#*=}"
            ;;
        esac
        case $i in
            era=*)
            era="${i#*=}"
            ;;
        esac
    done

    echo Found Proxy in: $X509_USER_PROXY
    if [ $isdata == "data" ]; then
        python WWG_postproc.py -d -y $year -e $era
    else
        python WWG_postproc.py -y $year -e $era
    fi

fi
