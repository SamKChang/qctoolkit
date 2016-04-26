import qctoolkit as qtk
import paramiko
import pexpect
import os
import shutil

def submit(inp_list, root, **remote_settings):
  necessary_list = [
    'ip',
    'submission_script',
  ]
  default_dict = {
    'username': None,
    'password': None,
    'flags': None,
    'n_cpu': 1,
  }
  for k, v in default_dict.iteritems():
    exec "%s = %s" % (k, v)
  if 'password' in remote_settings:
    password = remote_settings['password']
  if 'username' in remote_settings:
    username = remote_settings['username']
  if 'remote_path' not in remote_settings:
    remote_path = './%s' % root
  else:
    remote_path = remote_settings['remote_path']
  for s in necessary_list:
    if s not in remote_settings:
      qtk.exit('cluster setting:%s not defined' % s)
    else:
      exec "%s = '%s'" % (s, remote_settings[s])

  if type(inp_list) is not list:
    inp_list = [inp_list]
  program = inp_list[0].setting['program']

  if os.path.exists(root):
    qtk.warning("root directory %s exist, uploading existing folder"\
                % root)
  else:
    cwd = os.getcwd()
    os.makedirs(root)
    os.chdir(root)
    for inp in inp_list:
      inp.write(inp.molecule.name)
    os.chdir(cwd)

  paramiko_kwargs = {}
  if username:
    paramiko_kwargs['username'] = username
  if password:
    paramiko_kwargs['password'] = password
  ssh = paramiko.SSHClient()
  ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
  ssh.load_system_host_keys()
  ssh.connect(ip, **paramiko_kwargs)
  ssh_stdin, ssh_stdout, ssh_stderr = \
    ssh.exec_command('ls %s' % remote_path)
  sshout = ssh_stdout.read()
  if len(sshout) > 0:
    qtk.exit('remote path %s exists' % remote_path)

  ssh_newkey = 'Are you sure you want to continue connecting'
  patterns = [ssh_newkey, '[Pp]assword:', pexpect.EOF]
  if username:
    cmd = 'scp -qr %s %s@%s:%s' % (root, username, ip, remote_path)
  else:
    cmd = 'scp -qr %s %s:%s' % (root, ip, remote_path)

  qtk.report('submit', cmd)
  

  p = pexpect.spawn(cmd)
  i = p.expect(patterns)
  if i == 0:
    qtk.report('submit', 'adding %s to known_hosts' % ip)
    p.sendline('yes')
    i = p.expect(patterns)
  if i == 1:
    p.sendline(password)
    i = p.expect(patterns)
  if i == 2:
    if not p.before:
      qtk.report('submit', 'scp completed')
    else:
      qtk.warning('scp message: %s' % p.before)

  exe = qtk.setting.program_dict[program]
  remote_cmd = "%s %s %s %d '%s'" % \
    (submission_script, exe, remote_path, n_cpu, flags)
  ssh.exec_command(remote_cmd)
  qtk.report('submit', remote_cmd)
  ssh.close()
