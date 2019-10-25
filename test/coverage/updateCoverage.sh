#!/bin/bash
script_path=`realpath $0`
coverage_dir=`dirname ${script_path}`
test_dir=`dirname ${coverage_dir}`
package_dir=`dirname ${test_dir}`

# Move to coverage directory
cd ${coverage_dir}

# Get code coverage
coverage_rate=`./coverage.sh | grep -e "^TOTAL" | awk 'NF>1{print $NF}'`  # "100%"
coverage_rate=${coverage_rate::-1}  # Strip last character

# Get coverage color
coverage_color="blue"
if [ "${coverage_rate}" -lt 40 ]; then
    coverage_color="red"
elif [ ${coverage_rate} -lt 60 ]; then
    coverage_color="orange"
elif [ ${coverage_rate} -lt 70 ]; then
    coverage_color="yellow"
elif [ ${coverage_rate} -lt 80 ]; then
   coverage_color="yellowgreen"
elif [ ${coverage_rate} -lt 90 ]; then
    coverage_color="green"
elif [ ${coverage_rate} -le 100 ]; then
    coverage_color="brightgreen"
fi

# update README.md
old_badge_url="\(https:\/\/img.shields.io\/badge\/coverage-.+\%25-.+\)"
new_badge_url="\(https:\/\/img.shields.io\/badge\/coverage-${coverage_rate}\%25-${coverage_color}\)"
sed -i -E 's/'${old_badge_url}'/'${new_badge_url}'/g' ${package_dir}/README.md
