SYNASEARCH_HOME=/opt/synamatix/synasearch/synasearch-1.4.1-SNAPSHOT/

export PERL5LIB=$PERL5LIB:/opt/synamatix/synabase/synabase-v3.0.9/lib

# SynaBASE environment variailes
export SYNABASEROOT=/opt/synamatix/synabase/synabase-v3.0.9
export LD_LIBRARY_PATH=$SYNABASEROOT/lib
export LD_RUN_PATH=$SYNABASEROOT/lib:$LD_LIBRARY_PATH
export SXDBConfig=/opt/synamatix/synabase/sxdbconfig.xml
export SXDBROOT=/opt/synamatix/synabase/db

# Export PATH
export PATH=$PATH:$SYNABASEROOT/bin:$SYNASEARCH_HOME/bin

