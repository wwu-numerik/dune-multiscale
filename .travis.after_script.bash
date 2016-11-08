#!/bin/bash

if [[ $TRAVIS_JOB_NUMBER == *.1 ]] ; then
    git config --global hooks.clangformat ${CLANG_FORMAT}
    PYTHONPATH=${SUPERDIR}/scripts/python/ python3 -c "import travis_report as tp; tp.clang_format_status(\"${TRAVIS_BUILD_DIR}\")"
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make doc
    #./.travis/init_sshkey.sh $encrypted_6207618e9fa1_key $encrypted_6207618e9fa1_iv wwu-numerik.github.io
    #${SUPERDIR}/.travis/deploy_docs.sh ${MY_MODULE} ${DUNE_BUILD_DIR}
fi
