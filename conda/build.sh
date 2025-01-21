#!/bin/bash

mkdir -p ${PREFIX}/bin

target=${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}
mkdir -p ${target}

cp ${RECIPE_DIR}/*.py ${target}/
chmod +x ${target}/*.py
ln -s ${target}/pmlst.py ${PREFIX}/bin/pmlst.py
ln -s ${target}/pmlst.py ${PREFIX}/bin/pmlst

# copy script to download database
cp ${RECIPE_DIR}/download-db.sh ${target}/
chmod +x ${target}/download-db.sh
ln -s ${target}/download-db.sh ${PREFIX}/bin/download-db.sh

# set PMLST_DB variable on env activation
mkdir -p ${PREFIX}/etc/conda/activate.d ${PREFIX}/etc/conda/deactivate.d
cat <<EOF >> ${PREFIX}/etc/conda/activate.d/pmlst.sh
export PMLST_DB=${target}/pmlst_db/
EOF

cat <<EOF >> ${PREFIX}/etc/conda/deactivate.d/pmlst.sh
unset PMLST_DB
EOF
