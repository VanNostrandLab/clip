#!/usr/bin/env bash
set -e
unset PYTHONPATH

conda_install () {
  echo "    Installing $3 ...";
  conda install --channel="$1" --prefix="$2" --yes --quiet "$3" >/dev/null
  echo "    Successfully Installed $3.";
}

echo "Checking environment ..."
for package in conda wget git
do
  if command -v "${package}" >/dev/null; then
    echo "    ${package} ... installed";
  else
    echo "    Installing CLIP pipeline requires ${package} but it's not installed. Aborted!"; exit 1;
  fi
done
echo "Checking environment complete."

clip="$( cd .; pwd -P )"
eclip="${clip}/eclip"
merge_peak="${clip}/idr"
clip_venv="${clip}/venv"
merge_peak_venv="${merge_peak}/venv"

echo "Set up virtual environment for clip ...";

echo "    Installing Python (3.8) ...";
conda create --prefix="${clip_venv}" --yes --quiet python=3.8 >/dev/null;
echo "    Successful installed Python (3.8).";

for package in nodejs html5lib r-base cython
do
  conda_install "conda-forge" "${clip_venv}" "${package}"
done

for package in bedtools fastqc fastq-tools samtools=1.9 star=2.4.0j perl perl-app-cpanminus
do
  conda_install "bioconda" "${clip_venv}" "${package}"
done

echo "    Installing Perl packages ... "
"${clip_venv}/bin/cpanm" Statistics::Basic --quiet >/dev/null;
"${clip_venv}/bin/cpanm" Statistics::Distributions --quiet >/dev/null;
"${clip_venv}/bin/cpanm" Statistics::R --quiet >/dev/null;
echo "    Successful installed 3 Perl packages."
conda_install "mvdbeek" "${clip_venv}" "ucsc_tools";

echo "    Installing cwltool ..."
"${clip_venv}/bin/pip" install cwlref-runner --prefix="${clip_venv}" --quiet;
echo "    Successfully Installed cwltool."

for package in cutadapt umi_tools
do
  conda_install "bioconda" "${clip_venv}" "${package}"
done

echo "    Installing clipper ..."
git clone --quiet https://github.com/VanNostrandLab/clipper.git "${clip}/clipper";
"${clip_venv}/bin/pip" install "${clip}/clipper" --prefix="${clip_venv}" --quiet;
rm -rf "${clip}/clipper"
echo "    Successfully Installed clipper."

echo "    Installing eclipdemux ..."
git clone --quiet https://github.com/VanNostrandLab/eclipdemux.git "${clip}/eclipdemux";
"${clip_venv}/bin/pip" install "${clip}/eclipdemux" --prefix="${clip_venv}" --quiet;
rm -rf "${clip}/eclipdemux";
echo "    Successfully Installed eclipdemux."

echo "    Installing makebigwigfiles ..."
git clone --quiet https://github.com/YeoLab/makebigwigfiles.git "${clip}/makebigwigfiles";
"${clip_venv}/bin/pip" install "${clip}/makebigwigfiles" --prefix="${clip_venv}" --quiet;
rm -rf "${clip}/makebigwigfiles"
echo "    Successfully Installed makebigwigfiles."

echo "Setting up virtual environment for merge peak (IDR) pipeline in ${merge_peak} ..."
echo "    Installing Python (3.5.1) ...";
conda create --prefix="${merge_peak_venv}" --yes --quiet python=3.5.1 >/dev/null 1>&2;
echo "    Successful installed Python (3.5.1).";

#for package in cython numpy pandas scipy setuptools matplotlib
for package in cython numpy=1.11 pandas=0.20 scipy=0.18 setuptools=27.2 matplotlib=2.0
do
  conda_install "conda-forge" "${merge_peak_venv}" "${package}"
done

echo "    Installing IDR (2.0.2) ..."
wget --quiet https://github.com/nboley/idr/archive/2.0.2.zip;
unzip -qq 2.0.2.zip;
"${merge_peak_venv}/bin/pip" install "${clip}/idr-2.0.2" --prefix="${merge_peak_venv}" --quiet >/dev/null;
rm -rf idr-2.0.2 2.0.2.zip
echo "    Successfully Installed IDR (2.0.2)."
echo "Successfully set up virtual environment and installed merge peak (IDR) pipeline to ${merge_peak} ..."

echo "    Finalizing installation ..."
read -r -d '' clip_environment << EOF || true
#!/usr/bin/env bash

unset PYTHONPATH
export PATH="${clip_venv}/bin:\$PATH"
EOF
echo "${clip_environment}" > "${clip_venv}/clip_environment.sh"

read -r -d '' eclip_environment << EOF || true
#!/usr/bin/env bash

export PATH="${eclip}/bin:\$PATH"
export PATH="${eclip}/cwl:\$PATH"
export PATH="${eclip}/wf:\$PATH"
EOF
echo "${eclip_environment}" > "${clip_venv}/eclip_environment.sh"

read -r -d '' merge_peak_environment << EOF || true
#!/usr/bin/env bash

export PATH="${merge_peak_venv}/bin:\$PATH"
export PATH="${merge_peak}/bin:\$PATH"
export PATH="${merge_peak}/bin/perl:\$PATH"
export PATH="${merge_peak}/cwl:\$PATH"
export PATH="${merge_peak}/wf:\$PATH"
EOF
echo "${merge_peak_environment}" > "${clip_venv}/merge_peak_environment.sh"

read -r -d '' script << EOF || true
#!/usr/bin/env bash

source ${clip_venv}/clip_environment.sh
python ${clip_venv}/bin/clip.py \$@
EOF
echo "${script}" > "${clip}/clip"
chmod +x "${clip}/clip"
echo "    Successfully finalized installation."

clip_py="${clip}/source/clip.py"
bin_clip_py="${clip_venv}/bin/clip.py"
sed "s|ECLIP_ENVIRONMENT|${clip_venv}/eclip_environment.sh|" "${clip_py}" > "${bin_clip_py}"
sed "s|CLIP_ENVIRONMENT|${clip_venv}/clip_environment.sh|" "${clip_py}" > "${bin_clip_py}"
sed "s|ECLIP|${eclip}|g" "${clip_py}" > "${bin_clip_py}"
sed "s|MERGE_PEAK_ENVIRONMENT|${clip_venv}/merge_peak_environment.sh|" "${clip_py}" > "${bin_clip_py}"
sed "s|MERGE_PEAK|${merge_peak}|" "${clip_py}" > "${bin_clip_py}"

echo "Successfully installed and set up environment for CLIP pipeline."
echo "Run ${clip}/clip -h to see the usage."
