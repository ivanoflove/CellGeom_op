import subprocess

pro1 = subprocess.Popen(['python3', 'fluent.py'])
pro2 = subprocess.Popen(['python3', 'fluent1.py'])
pro3 = subprocess.Popen(['python3', 'fluent2.py'])
pro4 = subprocess.Popen(['python3', 'fluent3.py'])


pro1.wait()
pro2.wait()
pro3.wait()
pro4.wait()
