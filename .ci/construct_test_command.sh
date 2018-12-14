token=d9212b19-6bbd-4466-8f45-18414480f879
ci_env="bash <(curl -s https://codecov.io/env)"
ci_post_cov="bash <(curl -s https://codecov.io/bash) -t $token"
test_cmd="$ci_env && pytest /root/mantid_total_scattering/ && coverage run -m pytest /root/mantid_total_scattering/  && $ci_post_cov"

echo $test_cmd

