import qctoolkit as qtk
import sqlalchemy as q
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import or_
from datetime import datetime as dt
import os

Base = declarative_base()

class Entry(Base):
  __tablename__ = 'entries'

  id = q.Column(q.Integer, primary_key=True)
  datetime = q.Column(q.DateTime, nullable=False)
  content = q.Column(q.Text, nullable=False)
  comment = q.Column(q.Text)

  def __repr__(self):
    return "%s %s %s" % (self.datetime, self.content, self.comment)

class Logger(object):
  def __init__(self, path=':memory:', db_str = None, **kwargs):

    if not db_str:
      db_str = 'sqlite:///' + path

    self.engine = create_engine(db_str, **kwargs)
    self.name = db_str

    if os.path.exists(path):
      qtk.progress('DB', 'loading existing database: %s' % path)

    Base.metadata.create_all(self.engine)
    self.session = self.get_session(new=True)

  def push(self, content, comment=None):
    now = dt.now()
    if comment:
      entry = Entry(datetime=now, content=content, comment=comment)
    else:
      entry = Entry(datetime=now, content=content)

    self.session.add(entry)
    try:
      self.session.commit()
    except Exception as err:
      qtk.warning('can not commit, error: %s' % err)

  def list(self, filter_flag=''):
    filter_flag = r'%' + filter_flag + r'%'
    out = self.session.query(Entry).filter(
      or_(Entry.content.like(filter_flag), 
          Entry.comment.like(filter_flag)
      )
    ).all()
    return out

  def match(self, flag):
    out = self.session.query(Entry).filter(
      or_(Entry.content == flag,
          Entry.comment == flag)
    ).all()
    return out

  def commit(self):
    self.session.commit()


  def get_session(self, new=False):
    if not new:
      return self.session
    else:
      session = sessionmaker(bind=self.engine)
      return session()

