#
# Regular cron jobs for the sanepic package
#
0 4	* * *	root	[ -x /usr/bin/sanepic_maintenance ] && /usr/bin/sanepic_maintenance
