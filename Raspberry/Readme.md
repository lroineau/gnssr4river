# Raspberry Information
---
### Change the yaml config file

TODO, ask Roelof

### How to setup a service file

If you want your raspberry to automaticaly launch on boot a programm, you need to create a service file. This allows to simply plug your raspberry and directly execute a programm without having to open a terminal. This is very usefull for field survey.
Follow the step to create your service file. 

We first need to create a file for the script that will be used in the service file. If you are an advanced user in Linux, you can directly retrieve the file by cloning the git repository. Else, in the terminal, write 
```
vim name-of-your-file.py
``` 
and copy-paste the python script. Be aware of the path were your file be stored. Then make the file executable 

```
chmod a+x name-of-your-file.py
``` 

If you are not familiar with vim, here are the basics: 
* to go into instert mod, press ``i``
* to live insert mod, press ``esc``
* to leave the file and save it, press ``:wq``
* to leave the file without saving it, press ``:q``
    
Now to create the service, go into the right directory  
```
sudo /etc/systemd/system/name-of-your-service.service
```

Put the following into your service file and replace the lines with your [...] 
```
[Unit]
Description=your description
After=multi-user.target

[Service]
Type=simple
User=Your username
Restart=on-abort
ExecStart=/usr/bin/python path-to-your-file.py

[Install]
WantedBy=multi-user.target
```
Before starting the service, execute the following line. It reloads to take the change into account
```
sudo systemctl daemon-reload
```
You can then start the service, if you are in the right directory (``/etc/systemd/system/``)
```
sudo systemctl start your-service.service
```
You can see if your service file is properly running by executing
```
sudo systemctl status your-service.service
```
You can also of course stop the service
```
sudo systemctl stop your-service.service
```
Your service file is now running on your raspberry until the board is shutdown, and restart automaticaly when the board boots.

## Author

* **Lubin Roineau** [lubin_roineau](https://github.com/lroineau/)
