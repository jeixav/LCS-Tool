function publish_help
%PUBLISH_HELP Create LCS Toolbox help browser documentation

publish('animate_flow_help.m','evalCode',false);
publish('strain_lcs_script_help.m','evalCode',false);
publish('lcs_functions.m','evalCode',false);

s = warning('query','animate_flow:negative_delay');
warning('off','animate_flow:negative_delay')
publish('lcs_user_guide.m','evalCode',false);
warning(s)

publish('lcs_product_page.m','evalCode',false);

end
