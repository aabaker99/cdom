# R
sudo yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
sudo yum install -y R

# for biomaRt
sudo yum install -y libcurl-devel openssl-devel libxml2-devel

Rscript prep_instance.R
