#!/usr/bin/env bash
token=d9212b19-6bbd-4466-8f45-18414480f879
cd_repo="cd /root/mantid_total_scattering"
ci_env="bash <(curl -s https://codecov.io/env)"
run_pytest="mantidpython -m pytest"
run_main="mantidpython --classic ./bin/mantidtotalscattering --help"
run_coverage="coverage run -m pytest ."
ci_post_cov="bash <(curl -s https://codecov.io/bash) -t $token"

test_cmd="$cd_repo && $ci_env && $run_pytest && $run_main && $run_coverage && $ci_post_cov"

echo $test_cmd

